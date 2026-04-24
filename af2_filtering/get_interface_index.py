import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta import *
import argparse
import sys

logos_path = os.environ['LOGOS_PATH']
sys.path.append(os.path.join(logos_path))

from silent_tools import silent_tools

def get_interface_residues(silent_file, output_file):
    """
    Extract interface residues from a silent file and write to output.

    Args:
        silent_file: Path to the input silent file
        output_file: Path to the output file
    """
    pyrosetta.init('-in:file:silent_struct_type binary -beta_nov16 -holes:dalphaball /storage/d1/users/ruizhi/softwares/rosetta/rosetta.source.release-408/main/source/external/DAlpahBall/DAlphaBall.gcc')

    # Get all tags to process
    test = core.io.silent.SilentFileData(pyrosetta.rosetta.core.io.silent.SilentFileOptions())
    test.read_file(silent_file)

    tag_list = test.tags()

    with open(output_file, 'w') as out:
        for tag in tag_list:
            pose = pyrosetta.Pose()
            test.get_structure(tag).fill_pose(pose)
            pdb_info = pose.pdb_info()
            interface = pyrosetta.rosetta.protocols.interface.select_interface_residues(pose, 'A_B', 5)

            # Collect interface residues as a list
            interface_residues = []
            for i in range(len(interface)):
                j = i + 1
                if interface[j] == 1:
                    # Get true pdb index
                    interface_residues.append(pdb_info.pose2pdb(j))

            # Convert to comma-separated string
            residues_str = ','.join(interface_residues)

            # Output tag and residues on one line
            output_line = f"{tag}\t{residues_str}\n"
            out.write(output_line)
            #print(output_line.strip())

def main():
    parser = argparse.ArgumentParser(description='Extract interface residues from a Rosetta silent file')
    parser.add_argument('-i', '--input',
                        required=True,
                        help='Input silent file path')
    parser.add_argument('-o', '--output',
                        default='interface_residues.txt',
                        help='Output file path (default: interface_residues.txt)')

    args = parser.parse_args()

    get_interface_residues(args.input, args.output)

if __name__ == '__main__':
    main()
