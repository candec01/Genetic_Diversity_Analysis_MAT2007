#Import data set from pandas: 
import pandas as pd
import matplotlib.pyplot as plt
import math
#1. loading the data set
df = pd.read_csv(r"Data_set/coi_genetic_diversity_iucn.csv") 
df.info() #So we know how much the total data we have

#CLEANING
def cleaning_df(df): #all this steps will help us to trim and get only good data ready to use to plot, analyze and draw conclussion
    #missing data: 
    df.dropna(inplace=True) #with implace you are telling the program to rewrite those missing lines, over the original data frame. So they will be deleted.
    #in this data frame nothing was missing but, just so i know how to do it for future data sets.

    #irrelevant colums: 
    cols_categoricas = ["group", "class", "red_list_status", "category", "dataset"]
    cols_numericas = ["n", "length", "ntd_div_pi"]
    for col in cols_categoricas:
        print(f"column {col}: {df[col].nunique()} sublevels → {df[col].unique()}") #we will know how many subdivisions has every column, to delete if there is a column that has only one value.
   
    #if there is repeated rows:
    print(f"before: {df.shape}") 
    df.drop_duplicates(inplace=True)
    print(f"after: {df.shape}") #this is for me to know if there was any difference between the after and the before of the change.
    
    #chech if thre is typographic errors: 
    #just in case in the data sets there are spelling mistakes, columns that are the same but one is in capital letters and the other not... but same info... 
    fig,ax = plt.subplots(nrows=5, ncols=1, figsize=(30,30))
    fig.subplots_adjust(hspace=3)      
    for i, col in enumerate(cols_categoricas):
        counts = df[col].value_counts()
        ax[i].bar(counts.index, counts.values)
        ax[i].set_title(col)
        ax[i].tick_params(axis="x", rotation=45)
    plt.title("genetical data")
    plt.show() #seeing it helps you check for this mistakes, in this case non mistake was found, no not action was taken.

    #Create BINARY threat Status so I only have two options: "threatened" or "not-threatened":
    # Create an empty column
    df["threat_status"] = pd.NA 
    # Threatened species:
    # VU = Vulnerable
    # EN = Endangered
    # CR = Critically Endangered
    df.loc[df["red_list_status"].isin(['VU', 'EN', 'CR']),"threat_status"] = "Threatened" #.isin tells you for each row if the value of each columns it's in any of the categories you created. The program will select them (.lo) and gives them the 'lable' "Threatened" or whatever 
    # Not threatened species:
    # LC = Least Concern
    # NT = Near Threatened
    df.loc[df["red_list_status"].isin(['LC', 'NT']), "threat_status"] = "Not Threatened"
    # Check the result
    print("\nThreat status:")
    print(df['threat_status'].value_counts(dropna=False))
    
    # 7. CREATE DATASET FOR ANALYSIS
    # Remove species that cannot be classified into our two binary categories like DD, NE, EX, EW
    df_analysis = df.dropna(subset=['threat_status'])
    
    print("\nDataset for analysis:")
    print(df_analysis["threat_status"].value_counts())
    # Return the dataset ready for analysis
    return df_analysis

df_analysis = cleaning_df(df)

#Mann-Whitney U test: 
def mann_whitney_manual(group1, group2):
#First mix and cluster the data in to groups and rank them: 
    combined = pd.concat([group1, group2]).reset_index(drop=True)
    ranks = combined.rank()
    
    n1 = len(group1) #list them
    n2 = len(group2)

    #add up ranks from group 1: 
    R1 = ranks[:n1].sum()

    #Calculate U for each group using formula: 
    U1 = R1 - n1*(n1+1)/2
    U2 = n1*n2 - U1
    U = min(U1, U2)

    #Normal approximation for the p-value (valid because n1 and n2 are large)
    mean_U = n1 * n2 / 2
    std_U = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12) #fomulas 
    z = (U - mean_U) / std_U

    #p-value 
    p_value = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return U, z, p_value, n1, n2 

threatened = df_analysis.loc[df_analysis["threat_status"] == "Threatened", "ntd_div_pi"] #do the same thing as bf with the groups, but to tide it woth the nucleotide diversity.
not_threatened = df_analysis.loc[df_analysis["threat_status"] == "Not Threatened", "ntd_div_pi"] 
U, z, p_value, n1, n2  = mann_whitney_manual(threatened, not_threatened)
print(f"U statistic: {U}")
print(f"z-score: {z}")
print(f"p-value: {p_value}")
#thanks to this values, conclusions could be drawn.

#Boxplot: 
fig, ax = plt.subplots(figsize=(10,8))
data_to_plot = [not_threatened, threatened]
labels = ["Not Threatened", "Threatened"]
bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
colors = ["lightgreen", "orange"]  # green: Not Threatened; orange: Threatened. Making it pretty
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
font1 = {'family':'serif','color':'darkblue','size':20}
font2 = {'family':'serif','color':'darkred','size':15}
ax.set_title("Genetic Diversity (COI π) per Threat State", fontdict = font1)
ax.set_xlabel("Threat state", fontdict = font2 )
ax.set_ylabel("Nucleotide diversity (π)", fontdict = font2)
plt.show()

#Divide statistics:
print(df_analysis.groupby("threat_status")["ntd_div_pi"].describe())
#df_analysis.groupby('threat_status') this will divide in subgroups the value of the column (threat_status) creating groups : Threatened and Not-Threatened.  
#['ntd_div_pi']. This tells the program to only choose the data in this column. 
#.describe: for each group calculates the following information: count, mean, std, min, 25%, 50%, 75%, max.

N = n1 + n2  # 383 + 3974
effect_size_r = z / math.sqrt(N)
print(f"Effect size (r): {effect_size_r}")


#extra part: Figure out if there is a difference among the different classes (birds, mammals, insects...)
#first i have to fixed a Umbral, because maybe not every class has enought data for NT and T: 
#how many species per class?: 
n_per_class = df_analysis.groupby(["class", "threat_status"]).size().unstack(fill_value=0)
print(n_per_class)
#data obtained is higher than expected so there is no need, altough im going to run it to learn how to: 

#Now decide the umbral and cut: 
Min_Umbral = 20
#create a loop for each class: 
resultados = []

for clase in df_analysis["class"].unique():
    sub_df = df_analysis[df_analysis["class"] == clase] #this will create a loop for each subtimes (like reptiles, asd collect all the info in that sub_df)
    
    threatened_c = sub_df.loc[sub_df["threat_status"] == "Threatened", "ntd_div_pi"] #same here but with threatened...
    not_threatened_c = sub_df.loc[sub_df["threat_status"] == "Not Threatened", "ntd_div_pi"]
    
    # sieve: both groups have to have at leats Min_Umbral of species:
    if len(threatened_c) < Min_Umbral or len(not_threatened_c) < Min_Umbral:
        continue
    
    U, z, p_value, n1, n2 = mann_whitney_manual(threatened_c, not_threatened_c)
    N = len(threatened_c) + len(not_threatened_c)
    r = z / math.sqrt(N) #how big is the difference? 

    resultados.append({"class": clase,"n_threatened": len(threatened_c),"n_not_threatened": len(not_threatened_c),"median_threatened": threatened_c.median(),"median_not_threatened": not_threatened_c.median(), "U": U,"z": z,"p_value": p_value,"effect_size_r": r})

df_results = pd.DataFrame(resultados).sort_values("p_value")
print(df_results)

#FINAL PLOT: 
fig, ax = plt.subplots(figsize=(8,6))
colors = ["green" if p < 0.05 else "orange" for p in df_results["p_value"]]
font1 = {'family':'serif','color':'darkblue','size':20}
font2 = {'family':'serif','color':'darkred','size':15}
ax.barh(df_results['class'], df_results['effect_size_r'], color = colors)
ax.set_xlabel("Size Effect (r)", fontdict = font2)
ax.set_title("Genetic Diversity vs Risk of Extinction, per Taxonomic", fontdict = font1)
ax.axvline(0, color="black", linewidth=0.8)
for umbral, etiqueta in [(0.1, "small"), (0.3, "medium"), (0.5, "big")]:
    ax.axvline(-umbral, color="gray", linestyle=":", linewidth=0.7)
    ax.axvline(umbral, color="gray", linestyle=":", linewidth=0.7)
plt.tight_layout()
plt.show()











