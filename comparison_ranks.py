import pandas as pd
import matplotlib.pyplot as plt
from  matplotlib.pyplot import xlim
import seaborn as sns
import numpy as np
import json 



def df_organizer(csv1_original,csv2_mutant):
    """gets 2 csv files and orgenaize them into df"""
    original_df ,mutant_df= pd.read_csv(csv1_original),pd.read_csv(csv2_mutant)
    original_df,mutant_df=original_df[~original_df['peptide'].str.startswith('X')] ,mutant_df[~mutant_df['peptide'].str.startswith('X')]#removing redundat values
    return original_df,mutant_df


def dictionary_creator(df):
    """get df and return a dictionary that mapps key as the hla name and value as the
    list of values of the rank the dictionary values are in list for the usage of the graph creator later on"""
    new_dict={}
    for hla, rank in zip(df["HLA"],df["%Rank"]):
        if hla not in new_dict.keys():
            new_dict[hla]=[]
            new_dict[hla].append(rank)
        else:
            new_dict[hla].append(rank)
    return new_dict        


def graph_creator(binder_dictionary,original_dictionary):
    """gets the original dictionary and the binder dictionary and returns graphes
    that compare the rank to the refrence in the original"""
    for hla,rank in binder_dictionary.items():
        plt.clf() #clearing the plot in every itaration
        x=[x/x for x in rank] #create x values f 1 to each rank list , because x axis is categorial 
        plt.scatter(x,rank) #create graph    
        for i, txt in enumerate(rank): #gives anotation to x,y which is the rank
            plt.annotate(txt, (x[i], rank[i]))
        plt.xlabel(hla),plt.ylabel("% Rank") #giving titles to the axis
        plt.axhline(original_rank_values[hla]) #create base line from the refrence original peptide
        plt.savefig(hla +".png") #save the graph 
        
def full_deatailed_dictionary(df):
    """gets a df and return dictionary HLA as the key and the value is a tuple of the peptide name and rank"""  
    detailed_dictionary={}
    for hla,peptide,rank in zip(df["HLA"],df["peptide"],df["%Rank"]):#creating a dict filtering wb 
        if hla not in detailed_dictionary.keys():
            detailed_dictionary[hla]=[]
        detailed_dictionary[hla].append((peptide,rank))  
    return detailed_dictionary    

def orgenaize_the_chnages(peptide_changes,binder_dict,file_desired_name):
    """gets a dictionary with detailes about all the random changes and dictionaey with data about 
       the mutanats , their %rank and hla and return notpad detailed with the base change,the identity of the 
       peptide, and hla that it binds. we nedd to insert how we want to call to our file also returns a dictionary that map the changes and their influence of the all"""
    pep={}  
    detailed_file=open("C:/Users/Elinor/Desktop/תואר שני/{}.txt".format(file_desired_name),"w+")
    for hla in binder_dict.keys(): # for each hla in weak binders
        for peptide_rank_tup in binder_dict[hla]: #for each tuple of peptide rank in dict 
            if peptide_rank_tup[0] in peptide_changes.keys():#comparing peptides
               index_changes=peptide_changes[peptide_rank_tup[0]][1] #index of subsitutaion
               original_base=peptide_changes[peptide_rank_tup[0]][0] #the original base was
               base="".join((peptide_rank_tup[0][:index_changes],original_base,peptide_rank_tup[0][index_changes+1:]))
               what_happend="the {} change by {} the rank is {} the hla is {}"\
                .format(str(base),str(peptide_changes[peptide_rank_tup[0]]),peptide_rank_tup[1],hla)
               print(what_happend)
               detailed_file.write(what_happend+"\n")
    if peptide_changes[peptide_rank_tup[0]] not in pep.keys():
        pep[peptide_changes[peptide_rank_tup[0]]]=[]
        pep[peptide_changes[peptide_rank_tup[0]]].append(hla)
    pep[peptide_changes[peptide_rank_tup[0]]].append(hla)   
    return pep       
           
    
    
def finding_duplicants(full_detaailed_dictionary):
    """gets the full detailed dictionary (example shape- 'HLA-B8301': [('CPTINCERY', 0.7)) and return to list of 
       duplicants and the unike mutation , duplicants are the peptide that bind more than one allel""" 
    unike_peptides=[]
    dupli_peptides=[]
    for tup in full_detaailed_dictionary.values():
        for i in range(len(tup)):
            if tup[i][0] not in unike_peptides:
                unike_peptides.append(tup[i][0])
            else:
                 dupli_peptides.append(tup[i][0])    
    dupli_peptides= list(set(dupli_peptides)) 
    return [unike_peptides,dupli_peptides]
    
       
         
    

original_df=df_organizer("C:/Users/Elinor/Desktop/תואר שני/thursday.csv","C:/Users/Elinor/Desktop/תואר שני/the_one_peptide.csv")[1]
mutant_df=df_organizer("C:/Users/Elinor/Desktop/תואר שני/thursday.csv","C:/Users/Elinor/Desktop/תואר שני/the_one_peptide.csv")[0]    
original_rank_values=dictionary_creator(original_df)    
sb_df= mutant_df[mutant_df["%Rank"]<=0.5] #Creating df of strong binders
strong_binders=dictionary_creator(sb_df)
graph_creator(strong_binders,original_rank_values)

wb_df= mutant_df[(mutant_df["%Rank"] > 0.5) & (mutant_df["%Rank"] <= 2)]
wb_df_of_the_original= original_df[(original_df["%Rank"] > 0.5) & (original_df["%Rank"] <= 2)]
peptide_changes = np.load('C:/Users/Elinor/Desktop/תואר שני/random_changes_thursday.npy',allow_pickle='TRUE').item() #loading the text file of the mutants as dictionary

#for row,peptide,rank in zip(wb_df["HLA"],wb_df["peptide"],wb_df["%Rank"]):#creating a dict filtering wb 
#    weak_binders[row]=("wb",peptide,rank)  
weak_dict=dictionary_creator(wb_df)    #doing the same process to the weak binders
graph_creator(weak_dict,original_rank_values)
weak_binders=full_deatailed_dictionary(wb_df)
orgenaize_the_chnages(peptide_changes,weak_binders,"weak detailed")
strong_binders=full_deatailed_dictionary(sb_df)
orgenaize_the_chnages(peptide_changes,strong_binders,"strong detailed")



        

        
    



