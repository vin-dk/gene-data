import pandas as pd
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

num_genes = 50 

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

# PART 1 - MEAN CALC

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


class Gene:
    def __init__(self,gene_number, r_mean_score):
        self.gene = gene_number
        self.r_mean = r_mean_score
        
gene_list = []

class PerfectTricluster:
    def __init__(self, gene_index, exp_index, time_index):
        # gene, exp, time is (i,j,k) whose value is computed using the formula
        # this represents a single element
        self.gene = gene_index
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
                print(f"Calculation was {running_sum} + {tri_value_array[i].r}")
                running_sum += tri_value_array[i].r
                print(f"iter : {i} with j : {j}")
                print(f"Index is : {tri_value_array[i].gene},{tri_value_array[i].exp},{tri_value_array[i].time} at sum: {running_sum}")
                j += 1
                i += 1
            else: 
                r_mean = running_sum/12 
                gene = Gene(gene_num,r_mean)
                gene_list.append(gene)
                j = 0
                running_sum = 0
                gene_num += 1
                
                
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
    
PerfectTricluster.createGene(tricluster_object_elements)

for gene in gene_list:
    print(f"Gene: {gene.gene} and r mean: {gene.r_mean}")
    

# residuals calculated in (0,0,0)
# original minus perfect tricluster value calculated

# s for the tricluster is calculated : 1/IJK * sum(ALL RESIDUALS)^2, where each residual is squared -- THIS IS NOT IMPLEMENTED BECAUSE IT IS NOT IMPORTANT YET

# PART 3 GENE POOL 

# idea is to make blocks. I will allow user input for block selection. Follows this basic procedure
# Randomly generate number (with min >~5 max<~50% (??), that is how much goes into block 1. Take the mean of all r scores per gene, whose M(gene) governs that gene. 
# Select a threshold that will allow 75% (??) of values to remain in block. I don't know if this is worthwhile yet. Speak with Sudipta. Do note, some r scores negative, consider how to handle this.
# Any M(gene r scores)< threshold then gene is REMOVED from block, and added back to gene pool (these can both be arrays), the array storing blocks can be array of arrays
# Once this is done, block contains genes with r scores allowing for 75% to remain 
# this block is finalized, goes into block array
# Randomly generate new number, repeat the block process
# Once the amount of genes left falls below 15% (??), assign all remaining to the last block. This is to prevent blocks of only 2-3 genes. This can be changed.
# Process ends under 2 conditions : the final block has NO genes return to pool (skipping the step above), or the step above is executed, and the genes there are not returned to pool regardless of threshold

class Block: 
    def __init__(self, gene_amount):
        self.self = self.self

# PROBLEM : The machine seems INCREDIBLY PRECISE. Like remarkably so. Some of the sums are different numbers but lead to the EXACT same answer. Some are off by very very little. This makes sense if its the same machine

        













