import pandas as pd
import random
import numpy as np

# gene reference names as per sample dataset
gene_ref = []

# these values are recorded in "hour_condition" format
# a more general format, with expected entire matrix in one excel file will be created, this is for testing purposes / seeing the entire guts working
# we can pretend that these selected genes are members of a tricluster for the test
# the code in this file IS bloated and unneccessarily verbose. I want it to be as readable as possible and follow a path of logic

# Load Excel file
excel_file = r"C:\Users\13046\Desktop\data_set.xlsx"
df = pd.read_excel(excel_file)

# each column of excel sheet is an array, right now it is defined as "first 50 rows of column 0, first 50 rows of column 1..." 
# this follows the correct formatting for our purposes

algorithm_1 = False # by r mean
algorithm_2 = False # random block
algorithm_3 = True # chaos block 

num_genes = 4000

gene_ref = df.iloc[0:num_genes, 0].tolist()
zero_1 = df.iloc[0:num_genes, 1].tolist()
zero_2 = df.iloc[0:num_genes, 2].tolist()
zero_3 = df.iloc[0:num_genes, 3].tolist()
three_1 = df.iloc[0:num_genes, 4].tolist()
three_2 = df.iloc[0:num_genes, 5].tolist()
three_3 = df.iloc[0:num_genes, 6].tolist()
six_1 = df.iloc[0:num_genes, 7].tolist()
six_2 = df.iloc[0:num_genes, 8].tolist()
six_3 = df.iloc[0:num_genes, 9].tolist()
twelve_1 = df.iloc[0:num_genes, 10].tolist()
twelve_2 = df.iloc[0:num_genes, 11].tolist()
twelve_3 = df.iloc[0:num_genes, 12].tolist()

# list of relevant arrays
data_points = [zero_1, zero_2, zero_3, three_1, three_2, three_3, six_1, six_2, six_3, twelve_1, twelve_2, twelve_3]

# check for validity 
#for i, array in enumerate(data_points):
    #print(f"Array {i+1}: {array[:5]}")

exp_amount = 3

time_amount = 4

element_count = 0

gene_looper = 0 #these vars are for use in the perfect tricluster, check that class to see how they work

exp_looper = 0

time_looper = 0 

exp_index = 0

time_index = 0 # these are used to determine correct arrays

# these are values recorded for each gene_ref above. so we can consider it in this way that [1,3,1] refers to [gene 1, condition 3, zero hour] or [i,j,k) and so on
# I will be closely following the method outlined in the email, hand calculations should be easy to check logic


for data_point in data_points: # amount of elements in whole matrix, this can be changed if it is supposed to be amount of elements present in the calculation 
    for i in range(len(data_point)):
        element_count = element_count + 1

num_iterations_gene = len(zero_1)  # arbitrary array picked to get length of arrays

num_iterations_exp = exp_amount 

num_iterations_time = time_amount

mean_objects_gene = [] # every element of this array follows the gene_ref indexing

mean_objects_exp = []

mean_objects_time = [] 

mean_objects_tri = [] 

class MeanObject:
    # simple class to hold mean objects
    def __init__(self, iteration):
        self.iteration = iteration
        self.value = None

    def assign_value_gene(self, data_points):
        # M(i) calc
        # format : M(i) = exp/time * sum [i[c]] where c iterates DOWN the arrays above. Behaves as [1,1,1] + [1,2,1] + [1,3,1] + [1,2,2] + ...
        assert len(set(map(len, data_points))) == 1 # all arrays are same size
        sum_val = 0
        for data_point in data_points:
            i = self.iteration
            # print(f"element is : {element} at iteration {i}") - this ensures that we are tracking elements correctly
            sum_val = data_point[i] + sum_val
            self.value = (1/(exp_amount*time_amount)) * sum_val
    
    def assign_value_exp(self, data_points):
        #M(j) calc
        assert len(set(map(len, data_points))) == 1 # all arrays are same size
        sum_val = 0
        j = self.iteration
        for i in range(j,len(data_points),3):
            sum_val += sum(data_points[i])
        #print(f"sum val is : {sum_val}")
        #print(f"element count is : {element_count}")
        #print(f"time amount is : {time_amount}")
        self.value = (1/((element_count)*(time_amount))) * sum_val      
            
        
    def assign_value_time(self, data_points):
        #M(k) calc
        assert len(set(map(len, data_points))) == 1 # all arrays are same size
        sum_val = 0 
        k = self.iteration
        j = self.iteration * 3 
        for i in range(j,j+3):
            # print (f'Accessing array : {data_points[i]}, in group {j}')
            sum_val += sum(data_points[i])
        self.value = (1/((element_count)*(time_amount))) * sum_val
    
    def assign_value_tri(self, data_points):
        #M(ijk)
        assert len(set(map(len,data_points))) == 1 
        sum_val = 0
        for data_point in data_points:
            sum_val = sum(data_point) + sum_val
        self.value = (1/((element_count)*(time_amount)*(exp_amount))) * sum_val
        
class GeneRes:
    def __init__(self,gene_number, r_mean_score):
        self.gene = gene_number
        self.r_mean = r_mean_score
        
class PerfectTricluster:
    def __init__(self, gene_index, exp_index, time_index):
        # gene, exp, time is (i,j,k) whose value is computed using the formula
        # this represents a single element
        self.gene = gene_index
        self.id = gene_ref[gene_index]
        self.exp = exp_index
        self.time = time_index
        self.tri_mean = mean_objects_tri[0].value
        self.perfect_value = None
        self.real_value = None
        self.array_index = None
        self.r = None
    
    def compute_perfect_elements(self):
        # (m i J K + m I j K + m I J k - 2m I J K )
        gene = mean_objects_gene[self.gene]
        exp = mean_objects_exp[self.exp]
        time = mean_objects_time[self.time]
        self.perfect_value = (gene.value + exp.value + time.value - (2*self.tri_mean))
    
    def compute_r(self):
        # compare calculated value against real
        # since we already have indexes accounted for, we can use those to access the correct arrays
        if self.time == 0:
            if self.exp == 0:
                self.array_index = 0
            elif self.exp == 1:
                self.array_index = 1
            elif self.exp == 2:
                self.array_index = 2
        elif self.time == 1:
            if self.exp == 0:
                self.array_index = 3
            elif self.exp == 1:
                self.array_index = 4
            elif self.exp == 2:
                self.array_index = 5
        elif self.time == 2:
            if self.exp == 0:
                self.array_index = 6 # this is incredibly ugly, but works
            elif self.exp == 1:
                self.array_index = 7
            elif self.exp == 2:
                self.array_index = 8
        elif self.time == 3:
            if self.exp == 0:
                self.array_index = 9
            elif self.exp == 1:
                self.array_index = 10
            elif self.exp == 2:
                self.array_index = 11
        
        target_array = data_points[self.array_index]
        
        self.real_value = target_array[self.gene]
        
        self.r = self.real_value - self.perfect_value
    
    @staticmethod
    def createGene(tri_value_array):
        # our objects are still stored in (i,j,k format. We need them to be a gene object with relevant info) 
        # every gene has 12 entries. So we can compute the r mean easily by using this fact
        j = 0
        running_sum = 0
        gene_num = 1
        i = 0
        while i < ((12*num_genes) + 1):
            print(i,j)
            if j < 12:
                if i == ((12*num_genes)):
                    break
                #print(f"Calculation was {running_sum} + {tri_value_array[i].r}")
                running_sum += tri_value_array[i].r
                #print(f"iter : {i} with j : {j}")
                #print(f"Index is : {tri_value_array[i].gene},{tri_value_array[i].exp},{tri_value_array[i].time} at sum: {running_sum}")
                j += 1
                i += 1
            else: 
                r_mean = running_sum/12 
                gene = GeneRes(gene_num,r_mean)
                gene_list.append(gene)
                j = 0
                running_sum = 0
                gene_num += 1
         
# PART 1 - MEAN CALC

# GENE MEAN (i)           
for i in range (num_iterations_gene):
    # create mean objects, each representing a gene. For now it simply enumerates 1,2,3... which is associated with corresponding gene. In actual code it would be well defined
    mean_obj = MeanObject(iteration = i)
    mean_obj.assign_value_gene(data_points)
    mean_objects_gene.append(mean_obj)

for i, mean_obj in enumerate(mean_objects_gene):
    # this should print out the mean value associated with each mean object (in format mean_i = M(gene_i))
    print(f"Mean_{i+1} gene value: {mean_obj.value}")
    

# EXP MEAN (j)
for i in range (num_iterations_exp):
    mean_obj = MeanObject(iteration = i) 
    mean_obj.assign_value_exp(data_points)
    mean_objects_exp.append(mean_obj)
    
for i, mean_obj in enumerate(mean_objects_exp):
    print(f"Mean_{i+1} exp value: {mean_obj.value}")
    

# TIME MEAN (k)
for i in range (num_iterations_time):
    mean_obj = MeanObject(iteration = i) 
    mean_obj.assign_value_time(data_points)
    mean_objects_time.append(mean_obj)
    
for i, mean_obj in enumerate(mean_objects_time):
    print(f"Mean_{i+1} time value: {mean_obj.value}")
    
# TRICLUSTER MEAN (ijk)
mean_obj = MeanObject(iteration = 0)
mean_obj.assign_value_tri(data_points)
mean_objects_tri.append(mean_obj)

for i, mean_obj in enumerate(mean_objects_tri):
    print(f"Mean_{i+1} tri value: {mean_obj.value}")
    

# PART 2 - PERFECT SHIFTING FACTORS AND RESIDUALS
# Idea is to compute the perfect shifting tricluster value, per element, using the formula provided. Because of how the values are stored, it should be relatively easy to compute in a sequential format

tricluster_object_elements = []

gene_list = []
                       
for gene_looper in range(num_genes): # please note that the index is one less than expected due to starting at 0 
    # operates as (0,0,0) , (0,0,1) , (0,0,2) .... (0,1,0), (0,1,1) ... (1,0,0) ... and so on 
    exp_looper = 0
    
    while exp_looper < 3:
        time_looper = 0
        
        while time_looper < 4:
            tri_object = PerfectTricluster(gene_looper, exp_looper, time_looper)
            tri_object.compute_perfect_elements()
            tri_object.compute_r()
            tricluster_object_elements.append(tri_object)
            time_looper += 1
            
        exp_looper += 1
        
        if exp_looper == 3:
            break
        
        time_looper == 0
        
    
#for tri in tricluster_object_elements: 
    #gene_print = tri.gene + 1
    #exp_print = tri.exp + 1
    #time_print = tri.time * 3
    
    #if time_print == 9:
        #time_print = 12

    #print (f"Gene: {gene_print}, Exp: {exp_print}, Time: {time_print} with perfect value: {tri.perfect_value} and residue score: {tri.r}")
    
if algorithm_1 :
    
    PerfectTricluster.createGene(tricluster_object_elements)

    for gene in gene_list:
        print(f"Gene: {gene.gene} and r mean: {gene.r_mean}")
    

# residuals calculated in (0,0,0)
# original minus perfect tricluster value calculated

# s for the tricluster is calculated : 1/IJK * sum(ALL RESIDUALS)^2, where each residual is squared -- THIS IS NOT IMPLEMENTED BECAUSE IT IS NOT IMPORTANT YET

# PART 3 GENE POOL 

# idea is to make blocks. Follows this basic procedure: 
# start with gene 1 at (i,j,k) (1,1,1), compute shifting values ai + bj + ck + s. Using this, we generate a range of possible values, where the original value is the lower bound for for the range. So it may
# be : {12.5 - 15.5} for example. Then, iterate through ALL (i,j,k), any in range are assigned to that block, any out of range are excluded. Repeat this process with the out of range values, starting at
# first {i,j,k} there. once the generated range is produced, do the process again, which is block 2. Each block is an object with associated (i,j,k) values. the last 50 or so genes that are unrepresented,
# can make up the last block
# for the range of possible values, I will consider the range to be between the original value and the absolute value of the residual. So we will compute perfect tricluster (i,j,k) then calc residual
# the range will then be original(i,j,k) - perfect(i,j,k) 

total_genes = (element_count/12)
total_exp = 3
total_time = 4
total_coverage = total_genes * total_exp * total_time 

observed_exp = []
observed_time = []
observed_gene = []  

class Block: 
    def __init__(self, block, block_num):
        self.block = block
        self.block_num = block_num
        self.size = len(block)
        self.coverage = None
        self.TQI = None
        self.gene_count = None
        self.exp_count = None
        self.time_count = None
              
        
    def calcRange(self): 
        # displays highest + lowest values in block
        highest = 0
        lowest = 100000000
        for element in self.block:
            value = element.real_value
            if value < lowest:
                lowest = value
            if value > highest:
                highest = value
        # print(f"Highest value is: {highest} and lowest value is {lowest} across {self.size} number of elements at block {self.block_num + 1}")
        
    def calcCoverage(self):
        # calculates the coverage of a single tricluster
        # here, we consider num of genes to be count of single gene index, eg (1,0,0) and (1,0,1) would result in 1, not 2
        # same procedure goes for exp and time 
        
        # here is the problem, coverage is not meant to be implemented like this, I am aware. However, it is presupposed that some genes would be excluded. What we COULD do is,
        # get rid of the last block, of unfitting elements. We could play around with a formula to determine exclusion range. This would allow coverage to not always = 1 when applied
        # if we pretend that we have 5000 test genes, and 4500 of them get clustered (500 dont fit in anywhere, with all exp/time observed), then that yields coverage: 90, which is 
        # more in line with the results I was seeing from their observations. This is simply a place holder for logic atm
        
        # different idea, I will fill in coverage for single tricluster, but do it on a running total, which will then be compared against max
        
        my_observed_exp = []
        my_observed_time = []
        my_observed_gene = []
        
        for element in self.block:
             # at the end of this we should have array of unique indexes

            if element.gene not in observed_gene:
                observed_gene.append(element.gene)
            
            if element.exp not in observed_exp:
                observed_exp.append(element.exp)
            
            if element.time not in observed_time:
                observed_time.append(element.time)
                
                # again, extremely ugly, but effective
            
            if element.gene not in my_observed_gene:
                my_observed_gene.append(element.gene)
            
            if element.exp not in my_observed_exp:
                my_observed_exp.append(element.exp)
            
            if element.time not in my_observed_time:
                my_observed_time.append(element.time)
                
                
        self.gene_count = len(my_observed_gene)

        self.exp_count = len(my_observed_exp)

        self.time_count = len(my_observed_time)

        
        # numerator = self.gene_count * self.exp_count * self.time_count
        
        # self.coverage = ((numerator/total_coverage) * 100)

    
    def calcTQI(self):
        # calculates the triclustering quality index of a single tricluster
        # r should be computed for each element already, so we calculate the mean squared residue first
        r_list = []
        
        for element in self.block:
            r_list.append(element.r)
        
        mean = (sum(r_list))/(len(r_list))
        
        mean_squared = mean ** 2
        
        # I take it that volume is referring to standard volume (x*y*z), so I will calculate it as (num_gene * num_exp * num_time)
        
        volume = self.gene_count * self.exp_count * self.time_count
        
        self.TQI = (mean_squared/volume)
        
    # def calcSDB(self):
        # calculates statistical difference from background (of entire tricluster set) 
        # this is a place holder for discussion tmr
        
        
            
# we already have an array instatiated : tricluster_object_elements that has objects that have real value, r value, and indexing associated with it. The first element should be (0,0,0), so we'll start
# there, in this way, we can see blocks as just collections of these objects, because all their information is already self contained

if algorithm_2:
    
    blocks_done = False  # all blocks assigned if true
    all_blocks = []
    left_over = [] 
    constant = 0  # this can be changed but 0 for now
    block_amount = 0
    stop_at = 100  # same as above
    block_elements = tricluster_object_elements  # original full (i,j,k) list

    while not blocks_done:
        
        block = []  # a singular block
        left_over = []  # elements not fitting into block

        block_start = block_elements[0]  # start up with first element
        lower_bound = block_start.real_value
        upper_bound = lower_bound + abs(block_start.r) + constant

        for tri in block_elements:
            if tri.real_value >= lower_bound and tri.real_value <= upper_bound: 
                block.append(tri) # in range
            else:
                left_over.append(tri) # not in range
                
        all_blocks.append(Block(block, block_amount)) # in range get assigned into block

        block_amount += 1

        block_elements = left_over
        
        block = []

        if len(block_elements) <= stop_at: # handle terminating conditions, in this case we DO NOT append the last block
            for instance in block_elements:
                block.append(instance)
            # all_blocks.append(Block(block, block_amount))
            blocks_done = True

        elif not block_elements:
            blocks_done = True

        
    
    # analysis portion 
    length = 0
    total_coverage = 0
    tqi_sum = 0
    tqi_count = 0
    avg_tqi = 0
    
    for block in all_blocks:
        # important to calc coverage before TQI
        length += block.size
        block.calcRange()
        block.calcCoverage()
        block.calcTQI()
        tqi_count += 1 
        tqi_sum += block.TQI
    
    avg_tqi = tqi_sum/tqi_count
    
    denom = (total_genes * total_exp * total_time)
    numer = (len(observed_gene) * len(observed_exp) * len(observed_time))
    
    total_coverage = (numer/denom) * 100
    
    print(f"All blocks have coverage : {total_coverage} with average TQI of: {avg_tqi}")
        
    print(length)
    
    
    
if algorithm_3: 
    # follows much of a similar strategy as the above method, where we init blocks and cluster them accordingly. Here are a couple issues : Should paramaters (constant, terminating block size, etc)
    # be randomized as well? Additionally, we have no definitive terminating condition. I think initially, what I will do is do the algorithm as proposed, and go off average TQI (SDB still needs
    # implementing), as coverage will not be incredibly useful with the last block included. 
    
    # How i am implementing now is this: cluster everything as above, set its parameter - some amount as initial terminating condition. Then we randomly select a value, randomly manipulate it (0-1) and 
    # store the resultant blocks (with one different value) as a its own set. Then we compare tqi of both results. If original is better, we keep, if new is better, we discard. Do the process again. The 
    # process is over when some limit (satisfiablity is reached).
    
    # I am going to store all blocks, and compare access them via index, a collection of blocks is only added to the array if the tqi is better than before
    
    # important note, generations (the list), does not contain an actual copy of generations. That can be implemented easily, but for now, since I am using ONLY assignments / references to objects,
    # generations essentially contains the exact same list in each index. It is a placeholder for this implementation 
    
    def logSigmoid(x): # for brevity
        return 1 / (1 + np.exp(-np.log(x)))     
    
    def computeWeirdC(t, t_max, k): # this computes the log sig function. 
        random_num = random.uniform(0,1)
        log_sig = logSigmoid((0.5 * t_max - t)/(k))
        result = random_num * log_sig
        return result
    
    def computeNewX(x_old, weird_c): # this computes new element as paper describes
        random_num = random.uniform(0,1)
        x_new = (weird_c * random_num) + x_old
        return x_new
        
    
    done = False # determines termination of loop
    loop_count = 1 # how many time loop has ran
    best_tqi = 0 # current running best tqi
    original_tqi = 0 # this is good to have to compare against the best one we produced
    p2 = True  # future potential? I need to understand the algorithm more
    target = 1000000 # this is to compare tqi against, our goal is to get it lower
    generations = [] # this array holds all_blocks, which is an array that, at every index, contains an array representing a block. It is intended to hold every block structuring produced by algo
    max_gen = 500 # for now, this is the only stop condition. Arbitrarily selecting a threshold is meaningless. However, hopefully we can produce a lower tqi than the original
    k = 0.5 # I have no idea if this is a reasonable value for k
    
    
    last_element_value = 0  # this var is VERY important, as it is fundamental to how the loop works. it allows us to go back to the prior state if our idea was bad
    last_element = None # and save which object it was
    
    
    blocks_done = False  # all blocks assigned if true
    block_elements = tricluster_object_elements  # original full (i,j,k) list    
    
    
    while not done:
        
        all_blocks = []
        left_over = [] 
        constant = 0  # this can be changed but 0 for now
        block_amount = 0
        stop_at = 100  # same as above        
        
        while not blocks_done:
            
            total_genes = (element_count/12)
            total_exp = 3
            total_time = 4
            
            observed_exp_loop = []
            observed_time_loop = []
            observed_gene_loop = []              
            
            
            if loop_count > 1 and p2 == True: # we need to gen a new idea
                random_block = random.choice(all_blocks)
                random_element = random.choice(random_block.block)
                last_element = random_element
                # now we need to generate a new element to test the idea against 
                element_value = random_element.real_value
                last_element_value = element_value
                weird_c = computeWeirdC(loop_count, max_gen, k)
                new_x = computeNewX(element_value,weird_c)
                
                print(f"New x computed as : {new_x} compared against {last_element_value}")
                
                # now we have generated a new value for this randomly selected element. So we assign the element this value
                random_element.real_value = new_x
                # so now if we run the loop again, the value should shift ever so slighty. my only concern is that this only ever shifts values UP, because it cannot be negative
                

                
            block = []  # a singular block
            left_over = []  # elements not fitting into block

            block_start = block_elements[0]  # start up with first element
            lower_bound = block_start.real_value
            upper_bound = lower_bound + abs(block_start.r) + constant

            for tri in block_elements:
                
                if tri.gene not in observed_gene_loop:
                    observed_gene_loop.append(tri.gene)
                if tri.exp not in observed_exp_loop:
                    observed_exp_loop.append(tri.exp) # for calculating coverage
                if tri.time not in observed_time_loop:
                    observed_time_loop.append(tri.time)
                    
                if tri.real_value >= lower_bound and tri.real_value <= upper_bound: 
                    block.append(tri) # in range
                else:
                    left_over.append(tri) # not in range

            all_blocks.append(Block(block, block_amount)) # in range get assigned into block

            block_amount += 1

            block_elements = left_over

            block = []

            if len(block_elements) <= stop_at: # handle terminating conditions, in this case we DO NOT append the last block
                for instance in block_elements:
                    block.append(instance)
                all_blocks.append(Block(block, block_amount))
                blocks_done = True

            elif not block_elements:
                blocks_done = True


            # analysis portion 
            length = 0
            total_coverage = 0
            tqi_sum = 0
            tqi_count = 0
            avg_tqi = 0

            for block in all_blocks:
            # important to calc coverage before TQI
                length += block.size
                block.calcRange()
                block.calcCoverage()
                block.calcTQI()
                tqi_count += 1 
                tqi_sum += block.TQI    
            
            avg_tqi = tqi_sum/tqi_count
                
            denom = total_genes * total_exp * total_time
            numer = (len(observed_gene_loop) * len(observed_exp_loop) * len(observed_time_loop))
                
            total_coverage = (numer/denom) * 100
                
            print(f"All blocks have coverage : {total_coverage} with average TQI of: {avg_tqi}")
            
            if loop_count == 1:
                target = avg_tqi
                    
            if avg_tqi < target and loop_count > 1: # we add to our list if the idea is better
                generations.append(all_blocks)
                target = avg_tqi
                best_tqi = avg_tqi 
                print(f"New tqi accepted, {new_x} is better than {last_element_value}")
                
            elif avg_tqi >= target and loop_count > 1:
                last_element.real_value = last_element_value # reset the proper value
                print(f"New tqi rejected, {new_x} is worse. So the value is reset to it's original value of {last_element_value}")
                print(f"The target value was {last_element.perfect_value}")
                # here is the issue, if the "perfect value" is EVER lower, then the algorithm will never improve. Because it is strictly additive (on a always non-negative)
                
            else:
                original_tqi = avg_tqi
                
                
            
            loop_count += 1
            
            
            if loop_count <= max_gen:
                blocks_done = False
            elif loop_count == max_gen + 1: 
                done = True
            
            # now we have grouping as before, but we need to compare by choosing a random value
                        