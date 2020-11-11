import os
import optparse

parser = optparse.OptionParser()
parser.add_option("--min_ind", help="minimum index for the input output files", type=int)
parser.add_option("--max_ind", help="maximum index for the input output files", type=int)

(options, args) = parser.parse_args()


args = [str(i) for i in range(options.min_ind, options.max_ind + 1)]
for arg in args:
  os.system("python runRHAnalyzer.py --ind %s"%arg)

