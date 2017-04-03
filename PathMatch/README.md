# Run
python PathMatch.py example/query example/input example/corr result

# Taken from Pathmatch + Gaphmatch Paper

More information is available at http://faculty.cs.tamu.edu/shsze/pathmatch.

INPUT

The following files are needed:

1. A file that specifies the query path.  The names of vertices in the path
   are given in order.

   Example:

      Ste18p
      Ste4p
      Gpa1p
      Ste11p
      Ste5p
      Ste7p

   More examples are in ProteinInteractionNetwork/query and
   MetabolicNetwork/query.

2. A file that specifies the input graph in which related paths are sought.
   The line

      v1 v2 v3 v4

   defines edges (v1,v2), (v1,v3) and (v1,v4).  An undirected edge (v1,v2)
   should be included twice, once in v1 ... v2 ... and once in v2 ... v1 ...

   Example:

      ifc-2  dyb-1
      mig-15 ttx-1 dpy-14
      dyb-1  ifc-2
      dyp-14 mig-15
      ttx-1  mig-15

   A few pre-computed graphs from DIP (http://dip.doe-mbi.ucla.edu) and EcoCyc
   (http://ecocyc.org) are in ProteinInteractionNetwork/input and
   MetabolicNetwork/input.

3. A file specifying correspondences between vertices in the query path and
   the input graph.  Each correspondence is specified as

      v1 v2 corr

   where v1 is a vertex in the query path, v2 is a vertex in the input graph
   and corr is the similarity score.

   Example:

      Ste18p ifc-2  3e-04
      Ste11p mig-15 2e-03

   More examples are in ProteinInteractionNetwork/corr and
   MetabolicNetwork/corr.

USAGE

   python pathmatch query input corr output -g 1 -p 23.0 -n 3

positional arguments:
  query       file name containing query path
  input       file name containing input graph
  corr        file name containing correspondences between vertices
  output      output file name

optional arguments:
  -h, --help  show this help message and exit
  -g G        maximum number of mismatches or indels between two matches,
              default=1
  -p P        mismatch and indel penalty (>0), default=1
  -n N        number of output paths, default=3

OUTPUT

Each result is shown as a path alignment in which "--" between two vertices
denotes a match, "  " between two vertices denotes a mismatch and "-" denotes
an indel.  It is possible that some results are repeated since a path alignment
may be represented more than once.

Example:

Result01: score= 820.48
 Ste18p       -- ifc-2
   -             dyb-1
 Ste4p        -- unc-15
 Gpa1p           PIR:T24472
 Ste11p       -- mig-15
 Ste5p        -- ttx-1
 Ste7p        -- mig-15
   -             csn-5
 Fus3p        -- mpk-1
 Dig1p        -- GNB:17543358
 Ste12p       -- gei-4
 Mat1ap       -- PIR:T23537


