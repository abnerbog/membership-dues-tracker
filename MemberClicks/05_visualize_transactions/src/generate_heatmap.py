import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import gridspec, patches

def read_input_file(in_file):
    
    df = pd.read_csv(in_file,index_col=['Organization','Member Type','Initiation Year'])

    return df

def group_member_type(df):

    # graduate institutions
    df_gi = df.loc[df.index.get_level_values('Member Type') == 'University Member']
    # primarily undergraduate institutions
    df_pui = df.loc[df.index.get_level_values('Member Type') == 'Primarily Undergraduate Institution']
    # non-profit affiliate
    df_np = df.loc[df.index.get_level_values('Member Type') == 'Non-profit  Affiliate Members']
    # international affiliate
    df_intl = df.loc[df.index.get_level_values('Member Type') == 'International Affiliate Members']

    # create dictionary of dataframes
    dfs = {
        'University Member': df_gi,
        'Primarily Undergraduate Institution': df_pui,
        'Non-profit  Affiliate Members': df_np,
        'International Affiliate Members': df_intl
    }
  
    return dfs

def initialize_figure(dfs,N):

    # create figure and axes using gridspec
    plt.rcParams.update({'font.size': 4})
    fig = plt.figure(figsize=(4, 25))
    gs = gridspec.GridSpec(5, 1, height_ratios=[0.03, 
                                                len(dfs['University Member'])/N, 
                                                len(dfs['Primarily Undergraduate Institution'])/N, 
                                                len(dfs['Non-profit  Affiliate Members'])/N, 
                                                len(dfs['International Affiliate Members'])/N]
                                                )
    
    # initialize axes
    ax_legend = plt.subplot(gs[0])
    ax_gi = plt.subplot(gs[1])
    ax_pui = plt.subplot(gs[2])
    ax_np = plt.subplot(gs[3])
    ax_intl = plt.subplot(gs[4])

    axs = {
        'Legend':ax_legend,
        'University Member':ax_gi,
        'Primarily Undergraduate Institution':ax_pui,
        'Non-profit  Affiliate Members': ax_np,
        'International Affiliate Members': ax_intl
        }

    return fig,axs

def add_patch_based_on_value(row, col, value, ax, dues):
    # add light green rectangle if the dues are paid
    if value == dues:
        ax.add_patch(plt.Rectangle((col, row), 1, 1, color='#a3d1ac',alpha=1))
    # add dark green rectangle if there is an overpayment
    elif value > dues:
        ax.add_patch(plt.Rectangle((col, row), 1, 1, color='#5b9867',alpha=1))

def add_data(ax,df,dues):
    # create empty dataframe to create skeleton of heatmap
    df_empty = pd.DataFrame(0,index=df.index.get_level_values('Organization'),columns=df.columns)
    sns.heatmap(df_empty, ax=ax, cmap="Greys", cbar=False)
    # iterate through dataframe and add "cells" to heatmap where paymnets have been made
    for row in range(df.shape[0]):
        # grab the value in the initiation year index
        # will be empty string of length 1 if organization has not paid initiation fees 
        col_initiation = df.index.get_level_values('Initiation Year')[row]
        for col in range(df.shape[1]):
            # adds cell to heatmap if payment has been made
            add_patch_based_on_value(row, col, df.iat[row, col], ax, dues)
            # for initiation year add a red star 
            if df.columns[col] == col_initiation:
                ax.plot(col, row+0.5, 'r*', markersize=5)

def populate_heatmap(axs,dfs,current_year=2024):

    # note these values will change in future
    annual_payment_membertype = {
    'University Member': 200,
    'Primarily Undergraduate Institution': 70,
    'Non-profit  Affiliate Members': 100,
    'International Affiliate Members': 100
}

    # add data from each member type grouping to heatmap
    for member_type in dfs:
        add_data(axs[member_type],dfs[member_type],annual_payment_membertype[member_type])
        # create labels 
        axs[member_type].set_ylabel(member_type)
        # only set xticks and xlabels on the bottom grid
        if member_type == list(dfs.keys())[-1]:
            # set xlabels and xticks
            axs[member_type].set_xlabel('Membership Year')
            axs[member_type].set_xticks(axs[member_type].get_xticks()-0.5)
            axs[member_type].set_xticklabels(dfs[member_type].columns)
        else:
            axs[member_type].set_xlabel('')
            axs[member_type].set_xticks([]) 
        # add_vertical_line to current year for reference
        vline_pos = current_year - int(dfs[member_type].columns[0])
        line = axs[member_type].axvline(x=vline_pos, color='grey', linestyle='--', linewidth=1)
        line.set_dashes([5, 5])
        

def add_legend(ax_legend):

    # create indicator for initiation payment (red star)
    red_star_legend = mlines.Line2D([], [], color='white', marker='*', markerfacecolor='red', markersize=10, label='Paid (initiation)')
    # create indicator for paid/waived dues (light green cell)
    dues_2013_patch = patches.Patch(color='#a3d1ac', label='Paid or waived dues',alpha=1)
    # create indicator for prepay overpayment (dark green cell)
    dues_2024_patch = patches.Patch(color='#5b9867', label='Overpayment',alpha=1)

    # adding legend to legend axis
    ax_legend.legend(handles=[dues_2013_patch, red_star_legend, dues_2024_patch], loc='upper left', bbox_to_anchor=(0.1, 1), ncol=3)
    ax_legend.axis('off')  # hide the colorbar

def main(in_file,out_file):

    df = read_input_file(in_file)

    # define total number of institutions
    N = len(df)

    # group by member type on heatmap
    dfs = group_member_type(df)

    # initialize heatmap
    fig,axs = initialize_figure(dfs,N)

    # populate heatmap
    populate_heatmap(axs,dfs)

    # align y-labels in heatmap 
    fig.align_ylabels([axs[key] for key in axs if key != 'Legend'])

    # add legend
    add_legend(axs['Legend'])

    # save to file
    plt.tight_layout()
    plt.savefig(out_file)

if __name__ == '__main__':

    # inputs from snakefile
    in_filename = snakemake.input['in_filename']
    out_filename = snakemake.output['out_filename']

    main(in_filename,out_filename)
