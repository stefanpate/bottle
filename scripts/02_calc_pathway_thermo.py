from multiprocessing import Process, Queue
from src import thermo

# specify unique mongo connection string if using mongo
mongo_client_ID = "mongodb://ymc1840:yash_is_cool@minedatabase.ci.northwestern.edu:27017"

# specify expansion id
exp_ID = 'YC_00001'

# specify collection to write thermodynamic information to if using mongo
col_to_write = 'thermo'

num_processes = 5 # cores to use
mongo_client, pathway_numbers = thermo.get_relevant_pathways_frm_mongo(mongo_client_ID, exp_ID)

# initialize multiprocessing Queue
q = Queue()

# add all relevant pathway numbers to a shared multiprocessing Queue
for i, pathway_number in enumerate(pathway_numbers):
    q.put(pathway_number)

# close connection
mongo_client.close()

# define concentration lower and upper bounds (recommended by Dylan Brown of Lucks lab)
lb = 1e-6 # 1 micro-mol
ub = 500e-6 # 500 micro-mol

# start calculating pathway MDF values by multiprocessing
for i in range(num_processes):

    p = Process(target = thermo.calc_pathway_mdf, args=(q,
                                                        lb,
                                                        ub,
                                                        mongo_client_ID,
                                                        exp_ID,
                                                        col_to_write,)).start()




