import msprime
import numpy as np
import sys

def migration_example(chr, admix, t, file):
   
    gentime = 28
    #admix   = 0.0 #admixture proportion between groups

    #recombination map
    infile = "../../recomb_rates/HapmapII/genetic_map_GRCh37_chr" + str(chr) + ".txt.gz"
    recomb_map = msprime.RateMap.read_hapmap(infile)
    admix  = float(admix)

    #construct a demographic history
    #see https://tskit.dev/msprime/docs/stable/demography.html

    demography = msprime.Demography()
    
    #present-day pops
    Ne=5000
   
    pops = ["P1", "P2", "PX", "P3", "P4"]
    for i in range(0,4):
        demography.add_population(name='P' + str(i+1), initial_size=Ne)
    demography.add_population(name="PX", initial_size = Ne)

    demography.add_population(name="P12", initial_size=Ne)
    demography.add_population(name="P123", initial_size=Ne)
    demography.add_population(name="P1234", initial_size=Ne)
    print(pops)

    demography.add_population_parameters_change(time=0, population="P3", initial_size=1000);
    #define events (have to be ordered by time)
    demography.add_admixture(time = 50, derived = "PX", ancestral = ["P2", "P3"], proportions = [1-admix,admix])
    demography.add_population_parameters_change(time=50, population="P3", initial_size=Ne);
    demography.add_population_split(time=100, derived=["P1", "P2"], ancestral="P12")
    demography.add_population_split(time=int(t), derived=["P12", "P3"], ancestral="P123")
    demography.add_population_split(time=10000, derived=["P123", "P4"], ancestral="P1234")
 
    print(demography.debug())

    sam = []
    N = 10
    k = 0
    with open(file + "_tmp.fam", "w") as fp2:
        with open(file + ".poplabels", "w") as fp:
            fp.write("ID GROUP POP SEX\n")
            for p in pops:
              if p == "P3":
                sam.append(msprime.SampleSet(N, population = p))
                for i in range(0,N):
                  fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
                  fp2.write(p + " " + p + " " + p + " " + p + " 0 0\n")
                  k += 1
              elif p == "PX":
                sam.append(msprime.SampleSet(N, population = p))
                for i in range(0,N):
                  fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
                  fp2.write(p + " " + p + " " + p + " " + p + " 0 0\n")
                  k += 1
              elif p == "P1":
                sam.append(msprime.SampleSet(N, population = p))
                for i in range(0,N):
                  fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
                  fp2.write(p + " " + p + " " + p + " " + p + " 0 0\n")
                  k += 1
              elif p == "P4":
                sam.append(msprime.SampleSet(N, population = p))
                for i in range(0,N):
                  fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
                  fp2.write(p + " " + p + " " + p + " " + p + " 0 0\n")
                  k += 1
              else:
                sam.append(msprime.SampleSet(N, population = p))
                for i in range(0,N):
                  fp.write("tsk_" + str(k) + " " + p + " " + p + " NA\n")
                  fp2.write(p + " " + p + " " + p + " " + p + " 0 0\n")
                  k += 1
    fp.close()
    fp2.close()

    ts = msprime.sim_ancestry(
        samples = sam,
        ploidy = 2,
        demography = demography,
        recombination_rate = recomb_map#,
        #record_migrations = True
        )

    #add mutations
    mutated_ts = msprime.sim_mutations(ts, rate=1.25e-8, model = "binary")

    #write to file
    mutated_ts.dump(file + "_chr" + str(chr) + ".trees")

    n_dip_indv = int(ts.num_samples / 2)
    indv_names = [f"tsk_{str(i+1)}" for i in range(n_dip_indv)]
    with open(file + "_chr" + str(chr) + ".vcf", "w") as vcf_file:
        mutated_ts.write_vcf(vcf_file, individual_names=indv_names)

migration_example(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
