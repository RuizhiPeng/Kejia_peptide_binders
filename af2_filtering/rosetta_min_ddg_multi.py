import sys, os
import numpy as np
import pyrosetta
from collections import OrderedDict
import argparse
from multiprocessing import Pool, cpu_count
import tempfile
import shutil

logos_path = os.environ['LOGOS_PATH']
sys.path.append(os.path.join(logos_path))

try:
    from silent_tools import silent_tools
except ImportError:
    print("silent_tools not in path; adding directory from drhicks1")
    # The module isn't in the path, add the required directory and try again
    sys.path.append("/home/drhicks1/")
    try:
        from silent_tools import silent_tools
    except ImportError:
        # Handle the case where the module still can't be imported
        print("Failed to import silent_tools even after modifying sys.path")
        sys.exit(1)

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.std import ostringstream

from io import StringIO

pyrosetta.init( '-in:file:silent_struct_type binary  -ex1 -ex2 -beta_nov16 -holes:dalphaball  /storage/d1/users/ruizhi/softwares/rosetta/rosetta.source.release-408/main/source/external/DAlpahBall/DAlphaBall.gcc')

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

aa_1_N = {a:n for n,a in enumerate(alpha_1)}
aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}
aa_1_3 = {a:b for a,b in zip(alpha_1,alpha_3)}
aa_3_1 = {b:a for a,b in zip(alpha_1,alpha_3)}

fr_cart_fast_xml = f'''
<SCOREFXNS>
    <ScoreFunction name="sfxn" weights="beta_nov16" />
    <ScoreFunction name="sfxn_pack" weights="beta_nov16_cart" >
        <Reweight scoretype="coordinate_constraint" weight="1.0" />
        <Reweight scoretype="cart_bonded" weight="0.75" />
    </ScoreFunction>
    <ScoreFunction name="sfxn_min" weights="beta_nov16_cart" >
        <Reweight scoretype="coordinate_constraint" weight="0.1" />
        <Reweight scoretype="cart_bonded" weight="0.75" />
    </ScoreFunction>
</SCOREFXNS>

<RESIDUE_SELECTORS>
    <Chain name="chainA" chains="A"/>
    <Chain name="chainB" chains="B"/>
    <Neighborhood name="interface_chA" selector="chainB" distance="10.0" />
    <Neighborhood name="interface_chB" selector="chainA" distance="10.0" />
    <And name="AB_interface" selectors="interface_chA,interface_chB" />
    <Not name="Not_interface" selector="AB_interface" />

    <Slice name="chainA_last_res" indices="-1" selector="chainA" />
    <Slice name="chainB_first_res" indices="1" selector="chainB" />

    <And name="chainB_not_interface" selectors="Not_interface,chainB" />

    <And name="chainB_fixed" >
        <Or selectors="chainB_not_interface" />
    </And>
    <And name="chainB_not_fixed" selectors="chainB">
        <Not selector="chainB_fixed"/>
    </And>

    <True name="all" />

</RESIDUE_SELECTORS>

<TASKOPERATIONS>
    <IncludeCurrent name="current" />
    <ExtraRotamersGeneric name="ex1_ex2" ex1="1" ex2="1" />
    <LimitAromaChi2 name="limitchi2" chi2max="110" chi2min="70" include_trp="True" />

    <OperateOnResidueSubset name="restrict2repacking" selector="all">
        <RestrictToRepackingRLT/>
    </OperateOnResidueSubset>

    <OperateOnResidueSubset name="restrict_target_not_interface" selector="chainB_fixed">
        <PreventRepackingRLT/>
    </OperateOnResidueSubset>

</TASKOPERATIONS>

<MOVERS>
    <SwitchChainOrder name="chain1only" chain_order="1" />
    <SwitchChainOrder name="chain2only" chain_order="2" />
</MOVERS>
<SIMPLE_METRICS>
    <SapScoreMetric name="sap_score" score_selector="chainA" sap_calculate_selector="chainA" sasa_selector="chainA" />
</SIMPLE_METRICS>
<FILTERS>
    <ContactMolecularSurface name="contact_molecular_surface" distance_weight="0.5" target_selector="chainB" binder_selector="chainA" confidence="0" />
    <Ddg name="ddg_no_repack"  threshold="-10" jump="1" repeats="1" repack="0" confidence="0" scorefxn="sfxn" />
    <BuriedUnsatHbonds2 name="BUNS" jump_number="1"
        cutoff="20" generous_hbonds="true"
        sasa_burial_cutoff="0.01" AHD_cutoff="120"
        dist_cutoff="3.0" hxl_dist_cutoff="3.5"
        sulph_dist_cutoff="3.3" metal_dist_cutoff="2.7"
        scorefxn="sfxn"
        confidence="1.0" />
    <BuriedUnsatHbonds name="sbuns_interface" 
        residue_selector="AB_interface" 
        report_all_heavy_atom_unsats="true" 
        scorefxn="sfxn" 
        cutoff="4" 
        ignore_surface_res="false" 
        print_out_info_to_pdb="true" 
        use_ddG_style="true" 
        dalphaball_sasa="1" 
        burial_cutoff="0.01"
        burial_cutoff_apo="0.2"
        probe_radius="1.1" 
        max_hbond_energy="1.5"
        atomic_depth_selection="5.5" 
        atomic_depth_deeper_than="false" 
        confidence="0" />

    <BuriedUnsatHbonds name="Abuns" 
        residue_selector="all" 
        report_all_heavy_atom_unsats="true" 
        scorefxn="sfxn" 
        cutoff="4" 
        ignore_surface_res="false" 
        print_out_info_to_pdb="true" 
        use_ddG_style="true" 
        dalphaball_sasa="1" 
        probe_radius="1.1" 
        atomic_depth_selection="5.5" 
        atomic_depth_deeper_than="false" 
        confidence="0" />
</FILTERS>

<MOVERS>

    <ModifyVariantType name="remove_lower_terminus" remove_type="LOWER_TERMINUS_VARIANT,UPPER_TERMINUS_VARIANT" residue_selector="chainA_last_res" />
    <ModifyVariantType name="remove_upper_terminus" remove_type="UPPER_TERMINUS_VARIANT,LOWER_TERMINUS_VARIANT" residue_selector="chainB_first_res" />
    <AddChainBreak name="add_break" find_automatically="1" distance_cutoff="4"/>

    <PackRotamersMover name="cst_pack" scorefxn="sfxn_pack" task_operations="current,restrict2repacking,restrict_target_not_interface"/>
    <MinMover name="cart_min" max_iter="200" type="lbfgs_armijo_nonmonotone" tolerance="0.01" 
    cartesian="false" bondangle="true" bondlength="true" jump="1" bb="1" chi="1" scorefxn="sfxn_min" >
        <MoveMap name="MM"  >
            <Chain number="1" chi="true" bb="true" />
            <ResidueSelector selector="chainB_fixed" chi="false" bb="false" />
            <ResidueSelector selector="chainB_not_fixed" chi="true" bb="false" />
        </MoveMap>
    </MinMover>
</MOVERS>
'''

chainA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
chainB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")

objs = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(fr_cart_fast_xml)

cst_pack = objs.get_mover("cst_pack")
cart_min = objs.get_mover("cart_min")

termini1 = objs.get_mover("remove_lower_terminus")
termini2 = objs.get_mover("remove_upper_terminus")
add_break = objs.get_mover("add_break")

def add2scorefile(tag, scorefilename, write_header=False, score_dict=None):
        with open(scorefilename, "a") as f:
                add_to_score_file_open(tag, f, write_header, score_dict)

def get_final_dict(score_dict, string_dict) -> OrderedDict:
    '''
    Given dictionaries of numerical scores and a string scores, return a sorted dictionary
    of the scores, ready to be written to the scorefile.
    '''

    final_dict = OrderedDict()
    keys_score = [] if score_dict is None else list(score_dict)
    keys_string = [] if string_dict is None else list(string_dict)

    all_keys = keys_score + keys_string

    argsort = sorted(range(len(all_keys)), key=lambda x: all_keys[x])

    for idx in argsort:
        key = all_keys[idx]

        if ( idx < len(keys_score) ):
            final_dict[key] = "%8.3f"%(score_dict[key])
        else:
            final_dict[key] = string_dict[key]

    return final_dict

def add_to_score_file_open(tag, f, write_header=False, score_dict=None, string_dict=None):
        final_dict = get_final_dict( score_dict, string_dict )
        if ( write_header ):
                f.write("SCORE:  %s description\n"%(" ".join(final_dict.keys())))
        scores_string = " ".join(final_dict.values())
        f.write("SCORE:  %s    %s\n"%(scores_string, tag))

def add2silent( tag, pose, score_dict, sfd_out, output_filename="out.silent" ):
        # pose = pose_from_file( pdb )

        # pose = insert_chainbreaks( pose, binderlen )

        struct = sfd_out.create_SilentStructOP()
        struct.fill_struct( pose, tag )

        for scorename, value in score_dict.items():
            if ( isinstance(value, str) ):
                struct.add_string_value(scorename, value)
            else:
                struct.add_energy(scorename, value)

        sfd_out.add_structure( struct )
        sfd_out.write_silent_struct( struct, output_filename )

def record_checkpoint( tag_buffer, checkpoint_filename ):
        with open( checkpoint_filename, 'a' ) as f:
                for tag in tag_buffer:
                        f.write( tag )
                        f.write( '\n' )

def determine_finished_structs( checkpoint_filename ):
        done_set = set()
        if not os.path.isfile( checkpoint_filename ): return done_set

        with open( checkpoint_filename, 'r' ) as f:
                for line in f:
                        done_set.add( line.strip() )

        return done_set

def pose_from_silent(sfd_in, tag):
    pose = pyrosetta.Pose()
    sfd_in.get_structure(tag).fill_pose(pose)
    return pose

def get_filter_by_name(filtername):
    try:
        the_filter = objs.get_filter(filtername)
    except:
        the_filter = objs.get_simple_metric(filtername)
    # Get rid of stochastic filter
    if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
        the_filter = the_filter.subfilter()

    return the_filter

def filter_to_results(pose, filtername):
    this_filter = get_filter_by_name(filtername)
    if (isinstance(this_filter, pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter)):
        value = this_filter.compute(pose)
    else:
        value = this_filter.report_sm(pose)
    return value

def score_with_this_filter(pose, this_filter):
        return filter_to_results(pose, this_filter)

# harmonic near center. 0.98*max_val at radius
def my_topout_func(center, max_val, radius):
    limit = radius/2
    return pyrosetta.rosetta.core.scoring.func.TopOutFunc(max_val/limit**2, center, limit)


def generate_csts(pose, other_pose, subset, correct_rotamer_bonus, max_dev, visualize=True):

    ambiguous_pairs = {
        'E':[('OE1', 'OE2')],
        'D':[("OD1", 'OD2')],
        'F':[("CD1", 'CD2'),("CE1", 'CE2')],
        'Y':[("CD1", 'CD2'),("CE1", 'CE2')],
    }


    my_starts = []
    my_stops = []


    cst_set = pyrosetta.rosetta.utility.vector1_std_shared_ptr_const_core_scoring_constraints_Constraint_t()

    for seqpos in range(1, pose.size()+1):

        if ( not subset[seqpos] ):
            continue

        res = pose.residue(seqpos)
        other_res = other_pose.residue(seqpos)

        assert(res.name1() == other_res.name1())

        num_sc_atoms = res.nheavyatoms() - res.first_sidechain_atom() + 1

        this_weight = correct_rotamer_bonus / num_sc_atoms

        my_pairs = []
        if ( res.name1() in ambiguous_pairs ):
            my_pairs = ambiguous_pairs[res.name1()]

        # cs

        for i_sc_atom in range(num_sc_atoms):
            our_atom = res.first_sidechain_atom() + i_sc_atom

            if ( res.atom_type(our_atom).is_virtual() ):
                continue

            our_name = res.atom_name(our_atom)
            other_atom = other_res.atom_index(our_name)

            goal_xyz = other_res.xyz(other_atom)

            atom_id = pyrosetta.rosetta.core.id.AtomID( our_atom, seqpos )
            root = pyrosetta.rosetta.core.id.AtomID( 1, pose.size() )

            cst = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint( atom_id, root, goal_xyz, my_topout_func(0, this_weight, max_dev) )

            for my_pair in my_pairs:
                if ( our_name in my_pair ):
                    our_amb_atom_name = my_pair[1-my_pair.index(our_name)]

                    amb_cst = pyrosetta.rosetta.core.scoring.constraints.AmbiguousConstraint()
                    amb_cst.add_individual_constraint(cst)

                    amb_goal_xyz = other_res.xyz(our_amb_atom_name)

                    cst = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint( atom_id, root, amb_goal_xyz, my_topout_func(0, this_weight, max_dev) )

                    amb_cst.add_individual_constraint(cst)

                    cst = amb_cst

            cst_set.append(cst)

    return cst_set

def get_sap(pose):
    sap = pyrosetta.rosetta.core.pack.guidance_scoreterms.sap.SapScoreMetric()
    return sap.calculate(pose)

def pack_min_and_ddg(pose):
    correct_rotamer_bonus = 30
    max_dev = 2
    
    monomer_cst_set = generate_csts(pose, pose, chainA.apply(pose), correct_rotamer_bonus, max_dev, visualize=False)
    pose.add_constraints(monomer_cst_set)

    target_cst_set = generate_csts(pose, pose, chainB.apply(pose), correct_rotamer_bonus, max_dev, visualize=False)
    pose.add_constraints(target_cst_set)

    termini1.apply(pose)
    termini2.apply(pose)
    add_break.apply(pose)

    cst_pack.apply(pose)
    cart_min.apply(pose)
    cst_pack.apply(pose)
    cart_min.apply(pose)
    
    filters = ["ddg_no_repack", "contact_molecular_surface", "BUNS", "sbuns_interface", "Abuns"]
    filter_scores = {}
    for this_filter in filters:
        this_score = score_with_this_filter(pose, this_filter)
        filter_scores[this_filter] = this_score

    # Calculate SAP for chain A (binder) using the full complex pose
    filter_scores["sap_chainA"] = get_sap(pose.split_by_chain()[1])

    # Calculate SAP for all using the full complex pose
    filter_scores["sap_total"] = get_sap(pose)

    # Calculate ddg_per_sap, handling cases where sap_score is invalid
    if filter_scores["sap_chainA"] > 0 and filter_scores["sap_chainA"] < 900:
        filter_scores["ddg_per_sap"] = filter_scores["ddg_no_repack"] / filter_scores["sap_chainA"]
    else:
        filter_scores["ddg_per_sap"] = 999.0

    return pose, filter_scores

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
        # Initialize PyRosetta in this worker (if not already done)
        
        #pyrosetta.init('-in:file:silent_struct_type binary -beta_nov16 -mute all')

        # Reinitialize XML objects for this worker
        local_chainA = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("A")
        local_chainB = pyrosetta.rosetta.core.select.residue_selector.ChainSelector("B")
        local_objs = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(fr_cart_fast_xml)

        local_cst_pack = local_objs.get_mover("cst_pack")
        local_cart_min = local_objs.get_mover("cart_min")
        local_termini1 = local_objs.get_mover("remove_lower_terminus")
        local_termini2 = local_objs.get_mover("remove_upper_terminus")
        local_add_break = local_objs.get_mover("add_break")

        # Load structure from silent file
        sfd_in = pyrosetta.rosetta.core.io.silent.SilentFileData(
            pyrosetta.rosetta.core.io.silent.SilentFileOptions()
        )
        sfd_in.read_file(silent_file)

        pose = pyrosetta.Pose()
        sfd_in.get_structure(tag).fill_pose(pose)

        # Run pack_min_and_ddg with local objects
        correct_rotamer_bonus = 30
        max_dev = 2

        monomer_cst_set = generate_csts(pose, pose, local_chainA.apply(pose), correct_rotamer_bonus, max_dev, visualize=False)
        pose.add_constraints(monomer_cst_set)

        target_cst_set = generate_csts(pose, pose, local_chainB.apply(pose), correct_rotamer_bonus, max_dev, visualize=False)
        pose.add_constraints(target_cst_set)

        local_termini1.apply(pose)
        local_termini2.apply(pose)
        local_add_break.apply(pose)

        local_cst_pack.apply(pose)
        local_cart_min.apply(pose)
        local_cst_pack.apply(pose)
        local_cart_min.apply(pose)

        # Calculate filters
        filters = ["ddg_no_repack", "contact_molecular_surface", "BUNS", "sbuns_interface", "Abuns"]
        filter_scores = {}
        for this_filter in filters:
            local_filter = local_objs.get_filter(this_filter) if this_filter != "contact_molecular_surface" else local_objs.get_filter(this_filter)
            try:
                local_filter_obj = local_objs.get_filter(this_filter)
            except:
                local_filter_obj = local_objs.get_simple_metric(this_filter)

            if isinstance(local_filter_obj, pyrosetta.rosetta.protocols.filters.StochasticFilter):
                local_filter_obj = local_filter_obj.subfilter()

            if isinstance(local_filter_obj, pyrosetta.rosetta.protocols.simple_filters.ShapeComplementarityFilter):
                this_score = local_filter_obj.compute(pose)
            else:
                this_score = local_filter_obj.report_sm(pose)

            filter_scores[this_filter] = this_score

        # Calculate SAP for chain A (binder) using the full complex pose
        filter_scores["sap_chainA"] = get_sap(pose.split_by_chain()[1])

        # Calculate SAP for all using the full complex pose
        filter_scores["sap_total"] = get_sap(pose)
        filter_scores["delta_sap"] = filter_scores["sap_total"] - filter_scores["sap_chainA"]

        # Calculate ddg_per_sap, handling cases where sap_score is invalid
        if filter_scores["sap_chainA"] > 0 and filter_scores["sap_chainA"] < 900:
            filter_scores["ddg_per_sap"] = filter_scores["ddg_no_repack"] / filter_scores["sap_chainA"]
        else:
            filter_scores["ddg_per_sap"] = 999.0

        # Convert pose to PDB string for transfer
        oss = ostringstream()
        pose.dump_pdb(oss)
        pose_string = oss.str()

        print(f"Worker {worker_id}: Completed {tag}")
        return (f"{tag}_min", pose_string, filter_scores, True)

    except Exception as e:
        print(f"Worker {worker_id}: Error processing {tag}: {e}")
        import traceback
        traceback.print_exc()
        return (f"{tag}_min", None, None, False)

# If the script is executed directly, run a test or example function
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate ddG and SAP scores using multiprocessing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Use 8 cores
  python rosetta_min_ddg_multi.py input.silent -n 8

  # Use all available cores
  python rosetta_min_ddg_multi.py input.silent -n 0

  # Single core (default)
  python rosetta_min_ddg_multi.py input.silent
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

    if ncores == 1:
        # Single-threaded mode (original behavior)
        sfd_in = pyrosetta.rosetta.core.io.silent.SilentFileData(
            pyrosetta.rosetta.core.io.silent.SilentFileOptions()
        )
        sfd_in.read_file(silent)

        for idx, raw_jobname in enumerate(all_tags, 1):
            print(f"Processing {idx}/{len(all_tags)}: {raw_jobname}")
            pose = pose_from_silent(sfd_in, raw_jobname)
            pose, filter_scores = pack_min_and_ddg(pose)
            add2silent(f"{raw_jobname}_min", pose, filter_scores, sfd_out, f"{output_prefix}.silent")
            add2scorefile(f"{raw_jobname}_min", scorefilename, write_header=write_header, score_dict=filter_scores)
            write_header = False
            record_checkpoint([raw_jobname], checkpoint_filename)
    else:
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
