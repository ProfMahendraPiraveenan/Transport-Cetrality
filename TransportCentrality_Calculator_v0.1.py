import networkx as nx
import numpy as np
import math
import itertools
from itertools import count
import time
from time import sleep
import os
import collections
import tkinter
from tkinter import *
import tkinter.filedialog
import tkinter.simpledialog as simpledialog
import tkinter.messagebox
from decimal import Decimal
from statistics import mean
from scipy.stats import pearsonr


#main function that is called to compute Transport Centrality 
def runTC(input_file, output_file, beta):
    path, filename = os.path.split(input_file)
    print("\nINPUT: ", filename)
    
    edge_list,node_names = readFile(input_file)
    G = createNXGraph(edge_list,node_names)
    leaf_nodes = get_leafnodes(G)
    output_stat = output_file[:-4]+"_summary_stats.csv"
    f = open(output_stat, "w+")
    f.write("Beta Values, Average TC, Average BC, Ratio (AvgTC/AvgBC), Pearson Correlation(TC-BC)\n")
    f.close()
    for b0 in beta:
        b = round(b0,5)
        new_output_file = output_file[:-4]+"_beta"+str(b)+".csv"
        TC = usingNX_v2(G, b, leaf_nodes)
        
        norm_TC = normalizeTC_nx(TC,G)
        saveOutputToFile_nx(new_output_file, norm_TC, G, b, output_stat)

#rename raw nodes to ascending order starting from 0
def save_renamed_nodes(input_file, renamed_edges, original_edges):
    filename_w_ext = os.path.basename(input_file)
    filename_only, file_extension = os.path.splitext(filename_w_ext)
    path, filename = os.path.split(input_file)
    out_path = path+"/_RENAMED_nodes/"
    out_file = out_path+filename_only+"_renamedNodes.csv"

    try:
        os.makedirs(out_path)
    except FileExistsError:
        os.makedirs(out_path, exist_ok=True)
        pass

    f = open(out_file,"w+")
    for i in range(len(renamed_edges)):
        f.write("%d,%d,from,%d,%d\n" %(renamed_edges[i][0],renamed_edges[i][1],original_edges[i][0],original_edges[i][1]))


    before = np.array(original_edges).ravel()
    after = np.array(renamed_edges).ravel()
    

#get list of nodes with no children, or the leaf or terminal nodes
def get_leafnodes(G):
    leaf_dictionary = dict.fromkeys(G.nodes(),0)
    leaf_nodes = []
    is_directed_flag = nx.is_directed(G)
    for node in G.nodes():
        if is_directed_flag:        
            if(G.in_degree(node) == 0 or G.out_degree(node) == 0):
                leaf_nodes.append(node)
        else:
            if(len(G[node]) < 2):
                leaf_nodes.append(node)
    return leaf_nodes
    
#compute the exponential costs of paths given list of paths and the beta value
def computeTC_component(paths,  beta_value):
    cost_all = [len(path)-1 for path in paths]
    e_cost_all = [math.exp(-1*beta_value*cost) for cost in cost_all]
    paths_costs = sum(e_cost_all)
    return Decimal(paths_costs)

#transport centrality : finding the paths for each pair of nodes
def usingNX_v2(graph, beta_value, leaf_nodes):
    all_nodes = list(graph.nodes())
    size = graph.number_of_nodes()
    v_ratios = np.zeros((size, 1), dtype=Decimal)
    start = time.time()
    flatten = lambda l: [item for sublist in l for item in sublist]
    sum_cost = lambda paths: sum([len(path)-1 for path in paths])
    all_cost = lambda paths: [len(path)-1 for path in paths]
    e_cost = lambda all_cost: [math.exp(-1*beta_value*cost) for cost in all_cost]
    #sources = set(all_nodes) - set([v])
    for src in all_nodes:
        targets = set(all_nodes) - set([src])
        for target in targets:
            all_paths = list(nx.all_simple_paths(graph, src, target))
            
            if(len(all_paths) > 0):
                all_path_nodes = set(itertools.chain(*list(all_paths)))
                non_leaf_nodes  = set(all_path_nodes) - set(leaf_nodes)
                final_v = set(non_leaf_nodes) - set([src, target])
                for v in final_v:
                    v_paths = [path for path in all_paths if v in path]
                    
                    if (len(v_paths)> 0):
                        if(int(len(v_paths)) == int(len(all_paths))):
                            ratio = 1
                        else:

                            all_paths_costs = computeTC_component(all_paths, beta_value)
                            v_costs = computeTC_component(v_paths, beta_value)
                            ratio = v_costs / all_paths_costs
                        v_ratios[v] += ratio
                            
    print('\t*END of usingNX fxn w/execution time:', time.time() - start, '\n\n' )
   
    return v_ratios
     

#creates and returns a NetworkX graph object given list of edge pairs and node labels
def createNXGraph(edge_list, node_names):
    G = nx.Graph()
    G.add_edges_from(edge_list)
    G.remove_edges_from(G.selfloop_edges())
    print('=================================\n', nx.info(G),'\n=================================\n')
    temp = dict.fromkeys(G.nodes(), -1)
    nx.set_node_attributes(G, temp,'label')
    

    reference_names = dict()
    for i in range(len(node_names)):
        new_ = node_names[i][0]
        old_ = node_names[i][1]
        temp[new_] = old_
        reference_names[new_] = old_
    
    nx.set_node_attributes(G, temp,'label')

    return G
  
#reading files separated by commas (csv)
def readFile(filename):
    f = open (filename, 'r')
    lines = []
    for l in f:
        if(l.strip()):
            lines.append(l)
            
    l = [[num.strip() for num in line.split(',')] for line in lines]
    results = map(int, l)
    edge_list = np.zeros((len(l),2), dtype='int')
    node_dict = dict()
    node_names = []
    i = 0
    for r in range(len(l)):
      if(l[r]):
        try:
            src = int(l[r][0])
            dest = int(l[r][1])
            if node_dict.get(src, -1) == -1:   
                node_dict[src] = i
                node_names.append((i, src))
                i += 1
            if node_dict.get(dest, -1) == -1:
                node_dict[dest] = i
                node_names.append((i, dest))
                i += 1
                
            edge_list[r][0] = node_dict[src]
            edge_list[r][1] = node_dict[dest]
           
            
        except ValueError:
            pass
    return edge_list, node_names

#for writing output file separated by commas
def saveOutputToFile_nx(filename, norm_tc_list, graph, beta, summary_stat):
    node_list = list(graph.nodes())
    node_list.sort()
    f= open(filename,"w+")
    f_stat = open(summary_stat, 'a')
    
    bc = nx.betweenness_centrality(graph)
    bc_list = np.array(list(bc.values()))[:, np.newaxis]
    tc_list = list(norm_tc_list)
    
    r = pearsonr(bc_list.ravel(), norm_tc_list.ravel())

    #For the stats
    f.write("BETA VALUE,%.2f\n" %(beta))
    f.write("Average TC,%f\n" %(np.mean(norm_tc_list)))
    f.write("Average BC,%f\n" %(np.mean(bc_list))) # oh shit np.mean(list(np.array(bc)))
    

    temp_ratio = 0
    if mean(bc)==0:
        f.write("Ratio of Average (avg TC/ avg BC),%f\n" %(temp_ratio))
        f_stat.write("%.2f,%f,%f,%f,%f\n" %(beta, np.mean(norm_tc_list),np.mean(bc_list), temp_ratio, r[0]))
    else:
        temp_ratio = np.mean(norm_tc_list)/ np.mean(bc_list)
        f.write("Ratio of Average (avg TC/ avg BC),%f\n" %( np.mean(norm_tc_list)/ np.mean(bc_list) ))
        f_stat.write("%.2f,%f,%f,%f,%f\n" %(beta, np.mean(norm_tc_list),np.mean(bc_list),temp_ratio, r[0]))
    
    f.write("Pearson Correlation (TC-BC), %f\n\n" %(r[0]))
    #End of for the stats  
    f.write("Raw Node-ID(from inputdata), TC, BC , Ratio (TC/BC)")

    
    for node in node_list:
        if(bc[node]==0): #leaf node
            tc_bc_ratio = 0
            f.write("\n%d,%f,%f,%f" %(graph.node[node]['label'], norm_tc_list[node], bc[node], tc_bc_ratio))
        else:
            tc_bc_ratio = norm_tc_list[node] / bc[node]
            f.write("\n%d,%f,%f,%f" %(graph.node[node]['label'], norm_tc_list[node], bc[node], tc_bc_ratio))
    f.close()
    

#normalising the computed transport centrality given the calculated TC(v_ratios) and the Graph G
def normalizeTC_nx(v_ratios, G):
    all_nodes = list(G.nodes())  
    N = G.number_of_nodes()#len(all_nodes)
    norm = np.zeros((len(v_ratios),1), dtype='float')

    if nx.is_directed(G):
        print("G is directed")
        normalisation_factor = (N-1)
    else:
        print("G is undirected")
        normalisation_factor = (N-1)*(N-2)
        
    for i in range(len(v_ratios)):
         norm[i] = float(v_ratios[i])/ normalisation_factor

    return norm    

#loading nodes and grah
def loadData(input_data):  
    nodes = createNodes(input_data)
    graph = createGraph(input_data, nodes, 0)
    return nodes, graph


#function for initialising and createing the gui starts from here until the last line of code 
def gui():
    root = Tk()
    Title = root.title( "File Opener")
    label = ttk.Label(root, text ="I'm BATMAN!!!",foreground="red",font=("Helvetica", 16))
    label.pack()

    #Menu Bar

    menu = Menu(root)
    root.config(menu=menu)

    file = Menu(menu)

    file.add_command(label = 'Open', command = OpenFile)
    file.add_command(label = 'Exit', command = lambda:exit())

    menu.add_cascade(label = 'File', menu = file)



#function triggered for opening input file using single file processing 
def openFile():  
    batch_flag.set('False')
    current_dir = os.path.dirname(os.path.realpath(__file__))
    root.filename =  tkinter.filedialog.askopenfilename(initialdir = current_dir,title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
    input_text.set(root.filename)

#function triggered for opening input files using batch file processing 
def openBatchFiles():
    batch_flag.set('True')
    current_dir = os.path.dirname(os.path.realpath(__file__))
    root.directory =  tkinter.filedialog.askdirectory(initialdir = current_dir,title = "Select directory",)
    input_text.set(root.directory)

#function triggered for asking beta value 
def enterBeta():
    beta_val = simpledialog.askfloat("Enter beta value", "Beta", initialvalue=1)

    beta_text.set(beta_val)
    print("User beta value:", beta_val)

#function triggered for running TC using batch file processing 
def run_tc_batch():
    path = str(e1.get())
    user_beta = str(e2.get())
    beta_val_raw = round(float(user_beta), 5)
    beta_val = [beta_val_raw]

##    #To use with multiple beta values, uncomment code section below   
##    beta_val = np.arange(0.1,10.1,step=0.1) #from 1 to 10 EXPERIMENT 2.0

    new_path = path+"/_OUTPUT"
    files = []
    try:
        os.makedirs(new_path)
    except FileExistsError:
        # directory already exists
        os.makedirs(new_path, exist_ok=True)
        pass
    
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    
    for file in files:
        if file.endswith('.csv'):        
            output_text.set(new_path)
            
            curr_input = new_path+"/"+file
            out_temp = curr_input[0:-4]
            out_file = out_temp+"_OUTPUT_TC.csv"

            normalised_tc = runTC(path+"/"+file, out_file, beta_val)

#function triggered for running TC using single file processing 
def run_tc_single():
    
    user_beta = str(e2.get())
    beta_val_raw = round(float(user_beta), 5)
    beta_val = [beta_val_raw]
##    #To use with multiple beta values, uncomment code section below   
##    beta_val = np.arange(0.1,10.1,step=0.1)
    
##    tkinter.messagebox.showinfo("Output", str(batch_flag.get()))
    
    input_dir = os.path.split(root.filename)[0]
    input_file = os.path.split(root.filename)[1]
    new_path = input_dir+"/_OUTPUT/"
    
    
    try:
        os.makedirs(new_path)
    except FileExistsError:
        # directory already exists
        os.makedirs(new_path, exist_ok=True)
        pass
    
    out_file = input_file[0:-4]

    output_text.set(new_path+out_file+"_OUTPUT_TC.csv")
    normalised_tc = runTC(str(e1.get()), str(e3.get()), beta_val)
    
#function triggered when Run TC is clicked  
def run_tc():
    if (str(batch_flag.get()) == "False"):
        run_tc_single()
    else:
        run_tc_batch()
    
    
#function triggered when opening output using batch processing 
def openOutput_batch():
    
    path = str(e3.get())
    files = []
    message = 'Open output directory for output files'
    message+= "\n\nOUTPUT DIRECTORY: \n\n"+path
    message+= "\n\nOUTPUT FILES:\n\n"
    for (dirpath, dirnames, filenames) in os.walk(path):
        files.extend(filenames)
        break
    for file in files:
        message+=str(file)+"\n"
    tkinter.messagebox.showinfo("Output", message)


#function triggered when opening output using single file processing 
def openOutput():
    flag = str(batch_flag.get())
##    tkinter.messagebox.showinfo("Output", str(batch_flag.get()))
    if flag == 'True':
        openOutput_batch()
    elif flag=='False':
        tkinter.messagebox.showinfo("Output", str(e3.get()))
        filename = str(e3.get())
        try:
            with open(filename) as fp:
                message = fp.read()
        except IOError:
            message = 'No output available'
        tkinter.messagebox.showinfo("Output", message)
   
    
root = tkinter.Tk()
root.title("TC")

#variables global
input_text = tkinter.StringVar()
beta_text = tkinter.StringVar()
output_text = tkinter.StringVar()
batch_flag = tkinter.StringVar()


invisibleEntry = tkinter.Entry(root, text='Batch Flag', textvariable = batch_flag)

b0 = tkinter.Button(root, text='1. Batch Input Files', command=openBatchFiles )  
b1 = tkinter.Button(root, text='1. Choose Input File', command=openFile)  
e1 = tkinter.Entry(root, text='Input File', textvariable = input_text)
b2 = tkinter.Button(root, text='2. Enter Beta Value', command=enterBeta)  
e2 = tkinter.Entry(root, text='Enter Beta', textvariable = beta_text)


b4 = tkinter.Button(root, text='3. Compute TC', command=run_tc)  

l1 = tkinter.Label(root, text='TC results saved to output file below:')
e3 = tkinter.Entry(root, text='Output File', textvariable = output_text)
b3 = tkinter.Button(root, text='4. View Results', command=openOutput)

b0.pack()
b1.pack()
e1.pack(fill='x')
b2.pack(fill='x')
e2.pack(fill='x')
b4.pack(fill='x')

l1.pack(fill='x')
e3.pack(fill='x')
b3.pack(fill='x')

root.mainloop()

##def main():
##
##
##
##
##if __name__ == '__main__':
##    main()
    


    
