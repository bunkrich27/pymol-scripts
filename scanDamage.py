import sys, os
from pymol import cmd, stored, editing

@cmd.extend
def scanDamage( nucleosome, uvddb, cutoff=2, clashKeep=1, writeModels=False ):
  """
USAGE

  scanDamage nucleosome object,
             UVDDB probe object,
             distance in Ångstrom closer than which atoms are considered clashing,
             keep models with fewer than the specified number of clashes

  """
   
  chains = ['I','J']
  results=[]
  if writeModels:
    dirName =  None 
    dirName = ("%s_%s_models" % (nucleosome, uvddb) )
    if not os.path.exists(dirName):
      os.mkdir("%s" % (dirName) )
  
  for chain in chains: 

     chainLength = (len(cmd.get_model("polymer and %s and chain %s" % (nucleosome, chain)).get_residues()))

     for i in range(0, chainLength+1):
#        for i in range(130, 146):
        j=(i+1)
        
        if chain == 'I':
          comnplementaryChain = 'J'
        elif chain == 'J':
          matchedChain = 'I'


        # Handle lack of 5' phosphate on first base
        if i == 0: 
           cmd.select("moving", "%s and chain Z and resi 2 and backbone and not (n. P or n. OP1 or n. OP2)" % (uvddb) )
           cmd.select("target", "%s and chain %s and resi %d and backbone and not (n. P or n. OP1 or n. OP2)" % (nucleosome, chain, j) )
        elif i == 1: 
           cmd.select("moving", "%s and chain Z and resi 1-2 and backbone and not (n. P or n. OP1 or n. OP2)" % (uvddb) )
           cmd.select("target", "%s and chain %s and resi %d-%d and backbone and not (n. P or n. OP1 or n. OP2)" % (nucleosome, chain, i, j) )
        # Handle final base 
        elif i == chainLength:
           cmd.select("moving", "%s and chain Z and resi 1 and backbone" % (uvddb) )
           cmd.select("target", "%s and chain %s and resi %d and backbone" % (nucleosome, chain, i) )
        else:
           cmd.select("moving", "%s and chain Z and resi 1-2 and backbone" % (uvddb) )
           cmd.select("target", "%s and chain %s and resi %d-%d and backbone" % (nucleosome, chain, i, j) )

        #numAtomMoving=(cmd.count_atoms("moving"))
        #print ("\nNumber of moving atoms: %d" % (numAtomMoving) ) 
        #numAtomTarget=(cmd.count_atoms("target"))
        #print ("Number of target atoms: %d" % (numAtomTarget) ) 

        alignResult=[]

        alignResult = cmd.pair_fit("moving", "target")
        # alignResult = cmd.super("moving", "target")

        # print ("RMSD after refinement %.4f with %d atoms in alignhment after %d refinement cycles" % (alignResult[0],alignResult[1],alignResult[2]) )
        # print(" ".join('%s' % x for x in alignResult))

        clashes=("clashes")
        # Calculate clashes excluding DNA portion (chain Z resi 1-2) of UVDDB probe
        cmd.select("%s" % (clashes), "(%s and not chain Z) within %d of %s" % (uvddb, float(cutoff), nucleosome) )
        scoreAtoms = (cmd.count_atoms("%s" % (clashes) ))
        cmd.select("%s" % (clashes), "clashes byres clashes")
        scoreRes = (cmd.count_atoms("%s" % (clashes) + " and name CA"))

        # Write out models (nuclesome + current superimposed probe) for Rosetta scoring 
        if  writeModels :
          modelName = None
          modelName = ("%s_uvddb_%s_%d-%d" % (nucleosome, chain, i, j) )

          editing.copy_to("%s" % (modelName), "(%s or %s)" % (nucleosome, uvddb) )
          cmd.save("%s/%s.pdb" % (dirName,modelName), "(%s)" % (modelName) )
          cmd.delete("%s" % (modelName) )

        # Retain models with no. residue clashes less than or equal to 'clashKeep'
        if scoreRes < (int(clashKeep)+1):
           editing.copy_to("uvddb_%s_%d-%d" % (chain, i, j), "%s" % (uvddb) )
           cmd.group("ClashKeep %s)" % (clashKeep),"nonclashing + uvddb_%s_%d-%d" % (chain, i, j))
           cmd.order("uvddb_*", "yes")

        print ("DNA chain %s bases %d - %d :  UV-DDB clashes atoms %d residues %d " % (chain, i, j, scoreAtoms, scoreRes) )
        clashes=("%s,%d,%d,%d,%d" % (chain, i, j, scoreAtoms, scoreRes) )
        results.append((clashes)) 

  with open('%s_scanDamage.csv' % (nucleosome), 'w') as output:
    output.write("DNAchain,start,end,atomClashes,residueClashes\n")
    for i in results:
      output.write("%s\n" % i)
  print ("\nOutput saved to %s_scanDamage.csv\n" % (nucleosome) )
