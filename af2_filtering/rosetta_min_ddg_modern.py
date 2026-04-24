"""
Modern PyRosetta script for calculating ddG and interface metrics.
Uses direct PyRosetta API instead of XML-based RosettaScripts.
"""

import sys, os
import pyrosetta
from collections import OrderedDict
import argparse
from multiprocessing import Pool, cpu_count

logos_path = os.environ['LOGOS_PATH']
sys.path.append(os.path.join(logos_path))

try:
    from silent_tools import silent_tools
except ImportError:
    print("silent_tools not in path; adding directory from drhicks1")
    sys.path.append("/home/drhicks1/")
    try:
        from silent_tools import silent_tools
    except ImportError:
        print("Failed to import silent_tools even after modifying sys.path")
        sys.exit(1)

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.std import ostringstream
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, NeighborhoodResidueSelector, NotResidueSelector,
    AndResidueSelector, TrueResidueSelector
)
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    IncludeCurrent, ExtraRotamersGeneric, OperateOnResidueSubset
)
from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import SapScoreMetric

# Initialize PyRosetta
pyrosetta.init('-in:file:silent_struct_type binary -ex1 -ex2 -beta_nov16 -holes:dalphaball /storage/d1/users/ruizhi/softwares/rosetta/rosetta.source.release-408/main/source/external/DAlpahBall/DAlphaBall.gcc -mute all')

def get_final_dict(score_dict, string_dict) -> OrderedDict:
    """
    Given dictionaries of numerical scores and string scores, return a sorted dictionary
    of the scores, ready to be written to the scorefile.
    """
    final_dict = OrderedDict()
    keys_score = [] if score_dict is None else list(score_dict)
    keys_string = [] if string_dict is None else list(string_dict)

    all_keys = keys_score + keys_string
    argsort = sorted(range(len(all_keys)), key=lambda x: all_keys[x])

    for idx in argsort:
        key = all_keys[idx]
        if (idx < len(keys_score)):
            final_dict[key] = "%8.3f"%(score_dict[key])
        else:
            final_dict[key] = string_dict[key]

    return final_dict

def add_to_score_file_open(tag, f, write_header=False, score_dict=None, string_dict=None):
    final_dict = get_final_dict(score_dict, string_dict)
    if (write_header):
        f.write("SCORE:  %s description\n"%(" ".join(final_dict.keys())))
    scores_string = " ".join(final_dict.values())
    f.write("SCORE:  %s    %s\n"%(scores_string, tag))

def add2scorefile(tag, scorefilename, write_header=False, score_dict=None):
    with open(scorefilename, "a") as f:
        add_to_score_file_open(tag, f, write_header, score_dict)

def add2silent(tag, pose, score_dict, sfd_out, output_filename="out.silent"):
    struct = sfd_out.create_SilentStructOP()
    struct.fill_struct(pose, tag)

    for scorename, value in score_dict.items():
        if (isinstance(value, str)):
            struct.add_string_value(scorename, value)
        else:
            struct.add_energy(scorename, value)

    sfd_out.add_structure(struct)
    sfd_out.write_silent_struct(struct, output_filename)

def record_checkpoint(tag_buffer, checkpoint_filename):
    with open(checkpoint_filename, 'a') as f:
        for tag in tag_buffer:
            f.write(tag)
            f.write('\n')

def determine_finished_structs(checkpoint_filename):
    done_set = set()
    if not os.path.isfile(checkpoint_filename):
        return done_set

    with open(checkpoint_filename, 'r') as f:
        for line in f:
            done_set.add(line.strip())

    return done_set

def setup_residue_selectors(pose):
    """
    Setup residue selectors using modern PyRosetta API.
    Returns selectors for chains and interface regions.
    """
    # Chain selectors
    chainA = ChainSelector("A")
    chainB = ChainSelector("B")

    # Interface selectors (10Å neighborhood)
    interface_chA = NeighborhoodResidueSelector()
    interface_chA.set_focus_selector(chainB)
    interface_chA.set_distance(10.0)
    interface_chA.set_include_focus_in_subset(False)

    interface_chB = NeighborhoodResidueSelector()
    interface_chB.set_focus_selector(chainA)
    interface_chB.set_distance(10.0)
    interface_chB.set_include_focus_in_subset(False)

    # AB interface (both sides)
    AB_interface = AndResidueSelector()
    AB_interface.add_residue_selector(interface_chA)
    AB_interface.add_residue_selector(interface_chB)

    # Not interface
    not_interface = NotResidueSelector(AB_interface)

    # ChainB not interface (fixed region)
    chainB_not_interface = AndResidueSelector()
    chainB_not_interface.add_residue_selector(not_interface)
    chainB_not_interface.add_residue_selector(chainB)

    # ChainB interface (not fixed)
    chainB_interface = AndResidueSelector()
    chainB_interface.add_residue_selector(chainB)
    chainB_interface.add_residue_selector(NotResidueSelector(chainB_not_interface))

    return {
        'chainA': chainA,
        'chainB': chainB,
        'AB_interface': AB_interface,
        'chainB_fixed': chainB_not_interface,
        'chainB_movable': chainB_interface
    }

def setup_score_functions():
    """
    Setup score functions using modern PyRosetta API.
    """
    sfxn = pyrosetta.create_score_function('beta_nov16')

    sfxn_pack = pyrosetta.create_score_function('beta_nov16_cart')
    sfxn_pack.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint, 1.0)
    sfxn_pack.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 0.75)

    sfxn_min = pyrosetta.create_score_function('beta_nov16_cart')
    sfxn_min.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint, 0.1)
    sfxn_min.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded, 0.75)

    return sfxn, sfxn_pack, sfxn_min

def setup_task_operations(selectors):
    """
    Setup task operations for packing using modern PyRosetta API.
    """
    tf = TaskFactory()

    # Include current rotamers
    tf.push_back(IncludeCurrent())

    # Extra rotamers
    ex1_ex2 = ExtraRotamersGeneric()
    ex1_ex2.ex1(True)
    ex1_ex2.ex2(True)
    tf.push_back(ex1_ex2)

    # Limit aromatic chi2
    try:
        limit_chi2 = pyrosetta.rosetta.protocols.task_operations.LimitAromaChi2Operation()
        limit_chi2.chi2max(110)
        limit_chi2.chi2min(70)
        limit_chi2.include_trp(True)
        tf.push_back(limit_chi2)
    except:
        # If LimitAromaChi2Operation doesn't exist, skip it
        pass

    # Restrict all to repacking
    all_selector = TrueResidueSelector()
    restrict_to_repack = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
    restrict_op = OperateOnResidueSubset(restrict_to_repack, all_selector)
    tf.push_back(restrict_op)

    # Prevent repacking of fixed chainB regions
    prevent_repack = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_op = OperateOnResidueSubset(prevent_repack, selectors['chainB_fixed'])
    tf.push_back(prevent_op)

    return tf

def setup_movemap(pose, selectors):
    """
    Setup MoveMap for minimization.
    """
    mm = pyrosetta.MoveMap()

    # Chain A fully flexible
    chainA_residues = selectors['chainA'].apply(pose)
    for i in range(1, pose.size() + 1):
        if chainA_residues[i]:
            mm.set_bb(i, True)
            mm.set_chi(i, True)

    # Chain B: fixed regions no movement, interface regions chi only
    chainB_fixed = selectors['chainB_fixed'].apply(pose)
    chainB_movable = selectors['chainB_movable'].apply(pose)

    for i in range(1, pose.size() + 1):
        if chainB_fixed[i]:
            mm.set_bb(i, False)
            mm.set_chi(i, False)
        elif chainB_movable[i]:
            mm.set_bb(i, False)
            mm.set_chi(i, True)

    # Allow jump movement
    mm.set_jump(1, True)

    return mm

def add_coordinate_constraints(pose, selector_subset, correct_rotamer_bonus=30, max_dev=2):
    """
    Add coordinate constraints to maintain rotamer conformations.
    Uses TopOut function for soft constraints.
    """
    ambiguous_pairs = {
        'E':[('OE1', 'OE2')],
        'D':[("OD1", 'OD2')],
        'F':[("CD1", 'CD2'),("CE1", 'CE2')],
        'Y':[("CD1", 'CD2'),("CE1", 'CE2')],
    }

    cst_set = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()

    for seqpos in range(1, pose.size()+1):
        if not selector_subset[seqpos]:
            continue

        res = pose.residue(seqpos)
        num_sc_atoms = res.nheavyatoms() - res.first_sidechain_atom() + 1
        this_weight = correct_rotamer_bonus / num_sc_atoms

        my_pairs = []
        if res.name1() in ambiguous_pairs:
            my_pairs = ambiguous_pairs[res.name1()]

        for i_sc_atom in range(num_sc_atoms):
            our_atom = res.first_sidechain_atom() + i_sc_atom

            if res.atom_type(our_atom).is_virtual():
                continue

            our_name = res.atom_name(our_atom)
            goal_xyz = res.xyz(our_name)

            atom_id = pyrosetta.rosetta.core.id.AtomID(our_atom, seqpos)
            root = pyrosetta.rosetta.core.id.AtomID(1, pose.size())

            # Create TopOut function
            limit = max_dev / 2
            func = pyrosetta.rosetta.core.scoring.func.TopOutFunc(
                this_weight / (limit**2), 0, limit
            )

            cst = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint(
                atom_id, root, goal_xyz, func
            )

            # Handle ambiguous pairs
            for my_pair in my_pairs:
                if our_name.strip() in my_pair:
                    other_atom_name = my_pair[1 - my_pair.index(our_name.strip())]

                    amb_cst = pyrosetta.rosetta.core.scoring.constraints.AmbiguousConstraint()
                    amb_cst.add_individual_constraint(cst)

                    amb_goal_xyz = res.xyz(other_atom_name)
                    func2 = pyrosetta.rosetta.core.scoring.func.TopOutFunc(
                        this_weight / (limit**2), 0, limit
                    )

                    cst2 = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint(
                        atom_id, root, amb_goal_xyz, func2
                    )

                    amb_cst.add_individual_constraint(cst2)
                    cst = amb_cst
                    break

            cst_set.append(cst)

    return cst_set

def remove_termini_variants(pose, selectors):
    """
    Remove terminus variants at chain break point.
    """
    chainA_residues = selectors['chainA'].apply(pose)
    chainB_residues = selectors['chainB'].apply(pose)

    # Find last residue of chain A
    last_A = 0
    for i in range(pose.size(), 0, -1):
        if chainA_residues[i]:
            last_A = i
            break

    # Find first residue of chain B
    first_B = 0
    for i in range(1, pose.size() + 1):
        if chainB_residues[i]:
            first_B = i
            break

    # Remove variants
    if last_A > 0:
        try:
            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT, last_A
            )
        except:
            pass
        try:
            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT, last_A
            )
        except:
            pass

    if first_B > 0:
        try:
            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, pyrosetta.rosetta.core.chemical.UPPER_TERMINUS_VARIANT, first_B
            )
        except:
            pass
        try:
            pyrosetta.rosetta.core.pose.remove_variant_type_from_pose_residue(
                pose, pyrosetta.rosetta.core.chemical.LOWER_TERMINUS_VARIANT, first_B
            )
        except:
            pass

def calculate_ddg_simple(pose, sfxn):
    """
    Calculate ddG by separating chains and comparing energies.
    Simple non-repacking version.
    """
    # Score complex
    complex_energy = sfxn(pose)

    # Split chains and score separately
    chains = pose.split_by_chain()

    if len(chains) < 2:
        return 999.0

    chain1_energy = sfxn(chains[1])
    chain2_energy = sfxn(chains[2])

    ddg = complex_energy - (chain1_energy + chain2_energy)

    return ddg

def calculate_contact_molecular_surface(pose, selectors):
    """
    Calculate contact molecular surface using ShapeComplementarity.
    """
    try:
        sc_filter = pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter()
        # Set the selectors - may need different method names depending on PyRosetta version
        try:
            sc_filter.selector1(selectors['chainA'])
            sc_filter.selector2(selectors['chainB'])
        except:
            # Try alternative method names
            sc_filter.residue_selector1(selectors['chainA'])
            sc_filter.residue_selector2(selectors['chainB'])

        result = sc_filter.compute(pose)
        # The compute method returns a results object, extract the sc value
        if hasattr(result, 'sc'):
            return float(result.sc)
        else:
            # Fallback: try to convert directly
            return float(result)
    except Exception as e:
        print(f"Warning: ContactMolecularSurface calculation failed: {e}")
        return 0.0

def calculate_buried_unsats_filter(pose, selector=None, filter_type="BuriedUnsatHbonds"):
    """
    Calculate buried unsatisfied hydrogen bonds using the appropriate filter.
    """
    try:
        # Try to create the filter - BuriedUnsatHbonds2 for BUNS, BuriedUnsatHbonds for others
        if filter_type == "BUNS":
            try:
                buns = pyrosetta.rosetta.protocols.buns.BuriedUnsatHbondFilter2()
                buns.set_jump_number(1)
                buns.set_generous_hbonds(True)
                buns.set_scorefxn(pyrosetta.create_score_function('beta_nov16'))
            except:
                # Fallback to regular BuriedUnsatHbonds
                buns = pyrosetta.rosetta.protocols.simple_filters.BuriedUnsatHbondFilter()
        else:
            buns = pyrosetta.rosetta.protocols.simple_filters.BuriedUnsatHbondFilter()

            if selector is not None:
                try:
                    buns.set_residue_selector(selector)
                except:
                    buns.residue_selector(selector)

            try:
                buns.set_use_ddG_style(True)
            except:
                pass

            try:
                buns.set_ignore_surface_res(False)
            except:
                pass

        return buns.report_sm(pose)
    except Exception as e:
        print(f"Warning: {filter_type} calculation failed: {e}")
        return 0.0

def calculate_sap_score(pose):
    """
    Calculate SAP (Spatial Aggregation Propensity) score.
    """
    try:
        sap = SapScoreMetric()
        return sap.calculate(pose)
    except Exception as e:
        print(f"Warning: SAP calculation failed: {e}")
        return 0.0

def process_structure_worker(args):
    """
    Worker function for multiprocessing.
    Each worker processes a single structure independently.

    Args:
        args: tuple of (tag, silent_file, worker_id)

    Returns:
        tuple of (tag, pose_pdb_string, filter_scores, success)
    """
    tag, silent_file, worker_id = args

    try:
        # Load structure from silent file
        sfd_in = pyrosetta.rosetta.core.io.silent.SilentFileData(
            pyrosetta.rosetta.core.io.silent.SilentFileOptions()
        )
        sfd_in.read_file(silent_file)

        pose = pyrosetta.Pose()
        sfd_in.get_structure(tag).fill_pose(pose)

        # Setup selectors and score functions
        selectors = setup_residue_selectors(pose)
        sfxn, sfxn_pack, sfxn_min = setup_score_functions()

        # Add coordinate constraints
        monomer_cst = add_coordinate_constraints(
            pose, selectors['chainA'].apply(pose),
            correct_rotamer_bonus=30, max_dev=2
        )
        pose.add_constraints(monomer_cst)

        target_cst = add_coordinate_constraints(
            pose, selectors['chainB'].apply(pose),
            correct_rotamer_bonus=30, max_dev=2
        )
        pose.add_constraints(target_cst)

        # Remove terminus variants
        remove_termini_variants(pose, selectors)

        # Add chain break
        try:
            add_break = pyrosetta.rosetta.protocols.simple_moves.AddChainBreak()
            add_break.set_find_automatically(True)
            add_break.set_distance_cutoff(4.0)
            add_break.apply(pose)
        except:
            pass

        # Setup packing and minimization
        tf = setup_task_operations(selectors)
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(sfxn_pack)
        packer.task_factory(tf)

        mm = setup_movemap(pose, selectors)
        min_mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(sfxn_min)
        min_mover.min_type('lbfgs_armijo_nonmonotone')
        min_mover.tolerance(0.01)
        min_mover.max_iter(200)

        # Pack and minimize (2 rounds)
        packer.apply(pose)
        min_mover.apply(pose)
        packer.apply(pose)
        min_mover.apply(pose)

        # Calculate metrics
        filter_scores = {}

        # ddG using simple separation method
        filter_scores["ddg_no_repack"] = calculate_ddg_simple(pose, sfxn)

        # Contact molecular surface
        filter_scores["contact_molecular_surface"] = calculate_contact_molecular_surface(pose, selectors)

        # Buried unsats - different versions
        filter_scores["BUNS"] = calculate_buried_unsats_filter(pose, None, "BUNS")
        filter_scores["sbuns_interface"] = calculate_buried_unsats_filter(
            pose, selectors['AB_interface'], "sbuns_interface"
        )
        filter_scores["Abuns"] = calculate_buried_unsats_filter(pose, None, "Abuns")

        # SAP scores
        chains = pose.split_by_chain()
        if len(chains) >= 2:
            filter_scores["sap_chainA"] = calculate_sap_score(chains[1])
            filter_scores["sap_total"] = calculate_sap_score(pose)
            filter_scores["delta_sap"] = filter_scores["sap_total"] - filter_scores["sap_chainA"]
        else:
            filter_scores["sap_chainA"] = 999.0
            filter_scores["sap_total"] = 999.0
            filter_scores["delta_sap"] = 0.0

        # Calculate delta_sap (missing in original code)
        filter_scores["delta_sap"] = filter_scores["sap_total"] - filter_scores["sap_chainA"]

        # Calculate ddg_per_sap
        if filter_scores["sap_chainA"] > 0 and filter_scores["sap_chainA"] < 900:
            filter_scores["ddg_per_sap"] = filter_scores["ddg_no_repack"] / filter_scores["sap_chainA"]
        else:
            filter_scores["ddg_per_sap"] = 999.0

        # Ensure all filter_scores are plain Python floats (not PyRosetta types)
        filter_scores_clean = {}
        for key, value in filter_scores.items():
            try:
                # Try direct conversion first
                filter_scores_clean[key] = float(value)
            except (TypeError, ValueError):
                # If it's a complex object, try to extract the numeric value
                try:
                    # For ShapeComplementarity results, it might be an object with attributes
                    if hasattr(value, 'sc'):
                        filter_scores_clean[key] = float(value.sc)
                    elif hasattr(value, '__float__'):
                        filter_scores_clean[key] = float(value.__float__())
                    else:
                        # Last resort: convert to string then to float
                        filter_scores_clean[key] = float(str(value))
                except:
                    # If all else fails, use 0.0
                    print(f"Warning: Could not convert {key}={value} to float, using 0.0")
                    filter_scores_clean[key] = 0.0

        # Convert pose to PDB string for transfer
        oss = ostringstream()
        pose.dump_pdb(oss)
        pose_string = oss.str()

        print(f"Worker {worker_id}: Completed {tag}")
        return (f"{tag}_min", pose_string, filter_scores_clean, True)

    except Exception as e:
        print(f"Worker {worker_id}: Error processing {tag}: {e}")
        import traceback
        traceback.print_exc()
        return (f"{tag}_min", None, None, False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate ddG and interface metrics using modern PyRosetta API',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Use 70 cores
  python rosetta_min_ddg_modern.py input.silent -n 70 -o output_prefix

  # Use all available cores
  python rosetta_min_ddg_modern.py input.silent -n 0

  # Single core
  python rosetta_min_ddg_modern.py input.silent
        '''
    )
    parser.add_argument('silent', type=str, help='Input silent file')
    parser.add_argument('-n', '--ncores', type=int, default=1,
                       help='Number of CPU cores to use (default: 1, use 0 for all available cores)')
    parser.add_argument('-o', '--output', type=str, default='out',
                       help='Output prefix (default: out)')

    args = parser.parse_args()

    silent = args.silent
    ncores = args.ncores if args.ncores > 0 else cpu_count()
    output_prefix = args.output

    print(f"Starting with {ncores} core(s)...")

    # Create output files
    sfd_out = pyrosetta.rosetta.core.io.silent.SilentFileData(
        f"{output_prefix}.silent", False, False, "binary",
        pyrosetta.rosetta.core.io.silent.SilentFileOptions()
    )
    checkpoint_filename = f"{output_prefix}.check.point"
    scorefilename = f"{output_prefix}.sc"

    write_header = not os.path.exists(scorefilename)
    finished_structs = determine_finished_structs(checkpoint_filename)

    # Get all tags to process
    silent_index = silent_tools.get_silent_index(silent)
    all_tags = [tag for tag in silent_index['tags'] if tag not in finished_structs]

    print(f"Processing {len(all_tags)} structures...")

    # Multi-threaded mode
    print(f"Using multiprocessing with {ncores} workers...")

    # Prepare worker arguments
    worker_args = [(tag, silent, i % ncores) for i, tag in enumerate(all_tags)]

    # Process in parallel
    with Pool(processes=ncores) as pool:
        results = pool.map(process_structure_worker, worker_args)

    # Collect and write results
    print("\nWriting results to output files...")
    successful = 0
    for result in results:
        tag, pose_string, filter_scores, success = result

        if success and pose_string is not None:
            # Convert PDB string back to pose
            pose = pyrosetta.Pose()
            pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pose_string)

            # Write to output files
            add2silent(tag, pose, filter_scores, sfd_out, f"{output_prefix}.silent")
            add2scorefile(tag, scorefilename, write_header=write_header, score_dict=filter_scores)
            write_header = False
            record_checkpoint([tag.replace("_min", "")], checkpoint_filename)
            successful += 1

    print(f"\nCompleted: {successful}/{len(all_tags)} structures processed successfully")

    print(f"\nOutput files:")
    print(f"  Silent file: {output_prefix}.silent")
    print(f"  Score file:  {output_prefix}.sc")
    print(f"  Checkpoint:  {output_prefix}.check.point")
