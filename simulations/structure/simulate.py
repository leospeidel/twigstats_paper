import msprime
import numpy as np
import sys

def migration_example(chr):
   
    gentime = 28
    admix   = 0.25 #admixture proportion between groups

    #recombination map
    infile = "../recomb_rates/HapmapII/genetic_map_GRCh37_chr" + str(chr) + ".txt.gz"
    recomb_map = msprime.RateMap.read_hapmap(infile)

    #construct a demographic history
    #see https://tskit.dev/msprime/docs/stable/demography.html

    demography = msprime.Demography()
    
    #present-day pops
    pops = []
    for i in range(0,5):
      for j in range(0,5):
        demography.add_population(name='P' + str(i+1) + '_' + str(j+1), initial_size=500)
        pops.append('P' + str(i+1) + '_' + str(j+1))
    demography.add_population(name="anc", initial_size=500)
 
    print(pops)

    #define events (have to be ordered by time)
    for i in range(0,5):
      print(i)
      for j in range(0,5):
        if i >= 1 and i <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j+1), dest = 'P' + str(i)   + '_' + str(j+1), rate = 0.01)
        if i+2 >= 1 and i+2 <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j+1), dest = 'P' + str(i+2) + '_' + str(j+1), rate = 0.01)
        if j >= 1 and j <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j+1), dest = 'P' + str(i+1) + '_' + str(j),   rate = 0.01)
        if j+2 >= 1 and j+2 <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j+1), dest = 'P' + str(i+1) + '_' + str(j+2), rate = 0.01)
        if i >= 1 and i <= 5:
            demography.set_migration_rate( source = 'P' + str(i)   + '_' + str(j+1), dest = 'P' + str(i+1) + '_' + str(j+1), rate = 0.01)
        if i+2 >= 1 and i+2 <= 5:
            demography.set_migration_rate( source = 'P' + str(i+2) + '_' + str(j+1), dest = 'P' + str(i+1) + '_' + str(j+1), rate = 0.01)
        if j >= 1 and j <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j),   dest = 'P' + str(i+1) + '_' + str(j+1), rate = 0.01)
        if j+2 >= 1 and j+2 <= 5:
            demography.set_migration_rate( source = 'P' + str(i+1) + '_' + str(j+2), dest = 'P' + str(i+1) + '_' + str(j+1), rate = 0.01)

    demography.add_population_split(time=100, derived=pops, ancestral="anc")
 
    print(demography.debug())

    sam = []
    N = 10
    k = 0
    with open("msprime.poplabels", "w") as fp:
        fp.write("ID GROUP POP SEX\n")
        for p in pops:
          sam.append(msprime.SampleSet(N, population = p))
          for i in range(0,N):
            fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
            k += 1
    fp.close()

    ts = msprime.sim_ancestry(
        samples = sam,
        ploidy = 2,
        demography = demography,
        recombination_rate = recomb_map,
        record_migrations = True
        )

    #add mutations
    mutated_ts = msprime.sim_mutations(ts, rate=1.25e-8, model = "binary")

    #write to file
    mutated_ts.dump("msprime_chr" + str(chr) + ".trees")
    with open("msprime_chr" + str(chr) + ".vcf", "w") as vcf_file:
        mutated_ts.write_vcf(vcf_file)

migration_example(sys.argv[1])
