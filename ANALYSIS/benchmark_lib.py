### IMPORTS ###

from JSU_lib import *

import plotly.express as px
import plotly.graph_objects as go
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.proportion import proportion_confint

### FUNCTIONS ###

def get_best_preds_DIST(reference, predictor):
    """
    Gets best prediction for each reference binding site
    using the distance between centroids (minimum)
    """
    # Prepare data: Add suffixes to columns to prevent name clashes
    reference = reference.add_suffix('_ref')
    predictor = predictor.add_suffix('_pred')

    # Initialize the results list
    results = []

    # Iterate over each row in reference DataFrame
    for _, ref_row in reference.iterrows():
        # Filter predictor rows by 'rep_chain_pred' matching 'rep_chain_ref' of the reference row
        matched_predictor_rows = predictor[predictor['rep_chain_pred'] == ref_row['rep_chain_ref']] # subset rep_chain predictions
        
        if not matched_predictor_rows.empty:
            # Calculate Euclidean distances to the reference centre
            distances = matched_predictor_rows['centre_pred'].apply(lambda x: distance.euclidean(x, ref_row['centre_ref']))
            
            # Find the index of the minimum distance
            min_idx = distances.idxmin()
            min_distance = distances.min()
            
            # Append the row from predictor with closest centre
            closest_row = matched_predictor_rows.loc[min_idx]
            
            # Combine with the reference row and store
            combined_row = pd.concat([ref_row, closest_row])
            combined_row['min_distance'] = min_distance  # Store the minimum distance
        else: # NO predictions for a given structure
            nan_data = {col: np.nan for col in predictor.columns} # Create NaN row for predictors without matching rep_chain
            closest_row = pd.Series(nan_data)
            combined_row = pd.concat([ref_row, closest_row])
            combined_row['min_distance'] = np.nan  # No match found, so distance is NaN

        results.append(combined_row)

    return pd.DataFrame(results)

def get_best_preds_IREL(reference, predictor):
    """
    Gets best prediction for each reference binding site
    using the residue overlap (maximum)
    """
    # Prepare data: Add suffixes to columns to prevent name clashes
    reference = reference.add_suffix('_ref')
    predictor = predictor.add_suffix('_pred')

    # Initialize the results list
    results = []

    # Iterate over each row in reference DataFrame
    for _, ref_row in reference.iterrows():
        # Filter predictor rows by 'rep_chain_pred' matching 'rep_chain_ref' of the reference row
        matched_predictor_rows = predictor[predictor['rep_chain_pred'] == ref_row['rep_chain_ref']] # subset rep_chain predictions
        
        if not matched_predictor_rows.empty:

            def calculate_intersection(row):
                ref_set = set(ref_row['up_aas_ref'])
                pred_set = set(row['up_aas_pred'])
                return len(ref_set.intersection(pred_set)) / len(ref_set)
            
            # Calculate all relative intersections
            matched_predictor_rows['relative_intersection'] = matched_predictor_rows.apply(calculate_intersection, axis=1)
            
            # Filter rows with the highest relative intersection
            max_intersection = matched_predictor_rows['relative_intersection'].max()
            best_matches = matched_predictor_rows[matched_predictor_rows['relative_intersection'] == max_intersection]
            
            # If multiple matches have the same highest relative intersection, use Euclidean distance to choose the closest one
            if len(best_matches) > 1:
                best_matches['distance'] = best_matches['centre_pred'].apply(
                    lambda x: distance.euclidean(x, ref_row['centre_ref'])
                )
                best_match = best_matches.loc[best_matches['distance'].idxmin()]
            else:
                best_match = best_matches.iloc[0]
            
            # Store the best match row combined with reference row
            combined_row = pd.concat([ref_row, best_match])
        else: # NO predictions for a given structure
            nan_data = {col: np.nan for col in predictor.columns} # Create NaN row for predictors without matching rep_chain
            closest_row = pd.Series(nan_data)
            combined_row = pd.concat([ref_row, closest_row])
            combined_row['relative_intersection'] = np.nan  # No match found, so intersection is NaN

        results.append(combined_row)

    return pd.DataFrame(results)

def success_rate_DIST(df, dist_t_list, rank_t_list, total_predictions = None):

    results = []

    for origin_pred in df['origin_pred'].unique():
        subset = df[df['origin_pred'] == origin_pred] # subset of method predictions
        if total_predictions == None: # can pass absolute maximum (total observed pockets)
            total_predictions = len(subset) # otherwise, this will be the number of pockets in rows for which the method predicts
        print(f"Total predictions for {origin_pred} is {len(subset)}")

        for dist_t in dist_t_list: # for each distance threshold
            for rank_t in rank_t_list: # for each rank threshold
                condition = (subset['RANK_pred'] <= (subset['N_SITES_ref'] + rank_t)) & (subset['distance'] <= dist_t)
                correct_predictions = condition.sum()
                
                # Calculate the success rate
                success_rate = correct_predictions / total_predictions if total_predictions > 0 else 0
                
                # Calculate the standard error of the proportion
                se = (success_rate * (1 - success_rate) / total_predictions)**0.5 if total_predictions > 0 else 0
                
                # Calculate the confidence bounds
                lower_bound = max(0, success_rate - 1.96 * se)
                upper_bound = min(1, success_rate + 1.96 * se)
                
                results.append({
                    'origin_pred': origin_pred,
                    'rank_t': rank_t,
                    'dist_t': dist_t,
                    'success_rate': round(success_rate, 3),
                    'correct_predictions': correct_predictions,
                    'total': total_predictions,
                    'lower_se': round(lower_bound, 3),
                    'upper_se': round(upper_bound, 3)
                })

    return pd.DataFrame(results)

def success_rate_IREL(df, IRel_t_list, rank_t_list, total_predictions = None):
    results = []

    for origin_pred in df['origin_pred'].unique():
        subset = df[df['origin_pred'] == origin_pred]
        if total_predictions == None: # can pass absolute maximum (total observed pockets)
            total_predictions = len(subset) # otherwise, this will be the number of pockets in rows for which the method predicts
        print(f"Total predictions for {origin_pred} is {len(subset)}")

        for IRel_t in IRel_t_list: # for each IRel threshold
            for rank_t in rank_t_list: # for each rank threshold
                condition = (subset['RANK_pred'] <= (subset['N_SITES_ref'] + rank_t)) & (subset['relative_intersection'] >= IRel_t)
                correct_predictions = condition.sum()
                
                # Calculate the success rate
                success_rate = correct_predictions / total_predictions if total_predictions > 0 else 0
                
                # Calculate the standard error of the proportion
                se = (success_rate * (1 - success_rate) / total_predictions)**0.5 if total_predictions > 0 else 0
                
                # Calculate the confidence bounds
                lower_bound = max(0, success_rate - 1.96 * se)
                upper_bound = min(1, success_rate + 1.96 * se)
                
                results.append({
                    'origin_pred': origin_pred,
                    'rank_t': rank_t,
                    'IRel_t': IRel_t,
                    'success_rate': round(success_rate, 3),
                    'correct_predictions': correct_predictions,
                    'total': total_predictions,
                    'lower_se': round(lower_bound, 3),
                    'upper_se': round(upper_bound, 3)
                })

    return pd.DataFrame(results)

def get_volume_overlap(arr1: np.array, arr2: np.array, variant: str = 'jaccard') -> float:
    """
    Calculate the overlap between two sets of points represented by numpy arrays.

    Parameters:
    - arr1: np.array with shape (n, 3) representing points (x, y, z)
    - arr2: np.array with shape (n, 3) representing points (x, y, z)
    - variant: 'jaccard' for intersection over union, 'reference' for intersection over reference set

    Returns:
    - float: Overlap
    """
    if arr1.size == 0:
        return np.nan
    else:
        if arr2.size == 0:
            return 0.0
        else: # both have volume
            
            # Convert rows to tuples to represent points
            set1 = set(map(tuple, arr1))
            set2 = set(map(tuple, arr2))
            
            # Calculate intersection and union
            intersection = set1.intersection(set2)
            union = set1.union(set2)
            
            # Calculate Jaccard index based on the specified variant
            if variant == 'jaccard':
                jaccard = len(intersection) / len(union) if len(union) > 0 else 0.0
            elif variant == 'reference':
                jaccard = len(intersection) / len(set1) if len(set1) > 0 else 0.0
            else:
                raise ValueError("Invalid variant. Use 'jaccard' or 'reference'.")
            
        return round(jaccard, 4)

def compare_ALL_predictions(reference, predictor):
    """ Calculate closest rows from predictor for each row in reference using relative intersection and centroid proximity. """
    # Prepare data: Add suffixes to columns to prevent name clashes
    reference = reference.add_suffix('_ref')
    predictor = predictor.add_suffix('_pred')

    # Initialize the results list
    results = []

    # Iterate over each row in reference DataFrame
    for _, ref_row in reference.iterrows():
        # Filter predictor rows by 'rep_chain_pred' matching 'rep_chain_ref' of the reference row
        matched_predictor_rows = predictor[predictor['rep_chain_pred'] == ref_row['rep_chain_ref']]
        
        if not matched_predictor_rows.empty:

            def calculate_intersection(row):
                ref_set = set(ref_row['up_aas_ref'])
                pred_set = set(row['up_aas_pred'])
                return len(ref_set.intersection(pred_set)) / len(ref_set)
            
            # Calculate relative intersections and distances
            matched_predictor_rows['relative_intersection'] = matched_predictor_rows.apply(calculate_intersection, axis=1)
            matched_predictor_rows['distance'] = matched_predictor_rows['centre_pred'].apply(
                lambda x: distance.euclidean(x, ref_row['centre_ref'])
            )
            matched_predictor_rows['volume_overlap'] = matched_predictor_rows.apply(
                lambda row: get_volume_overlap(
                    volumes_dict.get(f"{ref_row['origin_ref']}", {}).get(f"{ref_row['rep_chain_ref']}_{ref_row['ID_ref']}", [None, np.array([])])[1],
                    volumes_dict.get(f"{row['origin_pred']}", {}).get(f"{row['rep_chain_pred']}_{row['ID_pred']}", [None, np.array([])])[1],
                    "reference"
                ), axis=1
            )
            
            # Combine reference row with each matched predictor row
            for _, pred_row in matched_predictor_rows.iterrows():
                combined_row = pd.concat([ref_row, pred_row])
                results.append(combined_row)
        else:
            # Create NaN row for predictors without matching rep_chain
            nan_data = {col: np.nan for col in predictor.columns}
            closest_row = pd.Series(nan_data)
            combined_row = pd.concat([ref_row, closest_row])
            combined_row['relative_intersection'] = np.nan  # No match found, so intersection is NaN
            combined_row['distance'] = np.nan  # No match found, so distance is NaN
            combined_row['volume_overlap'] = np.nan  # No match found, so volume overlap is NaN
            results.append(combined_row)

    return pd.DataFrame(results)

def normalise_score(group): # done to normalise IF-SitePred score
    if len(group) == 1:
        group['SCORE'] = 1
    else:
        max_val = group['n_points'].max()
        group['SCORE'] = group['n_points']/max_val
    return group

def plot_results_def(df, palette, markers_dict, fix_col, var_col, col_t, var_col_label, legend_order, ylabel = "Recall", xline = None, xticks=None, xlim = None, FSIZE = (10, 6), DPI = 100):
    # import matplotlib.pyplot as plt
    # import seaborn as sns

    plt.figure(figsize=FSIZE, dpi=DPI)

    df_col = df[df[fix_col] == col_t]

    sns.lineplot(data=df_col, x=var_col, y='success_rate', hue='origin_pred',
                 style='origin_pred', markers=markers_dict, palette=palette, 
                 markersize=10, linewidth=2.5, err_style=None, hue_order=legend_order, style_order=legend_order)

    for name, group in df_col.groupby('origin_pred'):
        plt.errorbar(group[var_col], group['success_rate'],
                     yerr=[group['success_rate'] - group['lower_se'], group['upper_se'] - group['success_rate']],
                     fmt='none', ecolor=palette[name], capsize=3)

    plt.xlabel(var_col_label)
    plt.ylabel(ylabel)
    leg = plt.legend(title='Method', title_fontsize='13', fontsize='11', loc='upper left', bbox_to_anchor=(1, 1), frameon=True)
    leg.get_frame().set_edgecolor('black')

    if xticks is not None:
        plt.xticks(xticks)

    if xlim != None:
        plt.xlim(xlim)

    if xline != None:
        plt.axvline(x = xline, linestyle = "--", color = "r", linewidth = 1)

    plt.show()

def jaccard_index(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0

def filter_redundant_pockets(df, threshold=0.5):
    # Initialize an empty list to store the filtered rows
    filtered_rows = []
    
    # Group by 'origin' and 'rep_chain'
    try:
        grouped = df.groupby(['origin', 'rep_chain'])
    except:
        grouped = df.groupby(['origin_pred', 'rep_chain_ref'])
    
    for (origin, rep_chain), group in grouped:
        group = group.reset_index(drop=True)
        
        # Create a mask to keep track of rows to keep
        keep_mask = np.ones(len(group), dtype=bool)
        
        for i in range(len(group)):
            if not keep_mask[i]:
                continue
            try:
                up_aas_i = set(group.loc[i, 'up_aas'])
            except:
                up_aas_i = set(group.loc[i, 'up_aas_pred'])
            try:
                score_i = group.loc[i, 'SCORE']
            except:
                score_i = group.loc[i, 'SCORE_pred']
            
            for j in range(i + 1, len(group)):
                try:
                    up_aas_j = set(group.loc[j, 'up_aas'])
                except:
                    up_aas_j = set(group.loc[j, 'up_aas_pred'])
                try:
                    score_j = group.loc[j, 'SCORE']
                except:
                    score_j = group.loc[j, 'SCORE_pred']
                
                jaccard = jaccard_index(up_aas_i, up_aas_j)
                
                if jaccard >= threshold:
                    if score_i >= score_j:
                        keep_mask[j] = False
                    else:
                        keep_mask[i] = False
                        break
        
        filtered_rows.extend(group[keep_mask].to_dict(orient='records'))
    
    return pd.DataFrame(filtered_rows)

def plot_cumulative_data(df, figsize=(6, 6), dpi=150, xlim=(0, 100), ylim=(0, 500), methods = "ALL", out = None):
    # Determine linestyle based on the presence of 'PRANK' in the 'origin_pred' column
    plt.figure(figsize=figsize, dpi=dpi)
    if methods != "ALL":
        df = df.query('origin in @methods')
    for method in methods:
        subset = df.query('origin == @method')
        #print(method, linestyles_dict2[method])
        plt.plot(subset['cumsum_FP'], subset['cumsum_TP'], label=method, color=palette[method], lw=1.5, linestyle = linestyles_dict2[method])

    plt.xlabel('# FP')
    plt.ylabel('# TP')
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    #plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.legend().set_visible(False)
    if out != None:
        plt.savefig(out)
    
    plt.show()

def plot_precision_data(df, figsize=(6, 6), dpi=150, xlim=(-1, 500), ylim=(0, 101), methods = "ALL", out = None):
    plt.figure(figsize=figsize, dpi=dpi)
    if methods != "ALL":
        df = df.query('origin in @methods')
    for method in methods:
        subset = df.query('origin == @method')
        plt.plot(
            subset['total'], subset['precision'], label=method, color=palette[method], lw=1.5,
            linestyle = linestyles_dict2[method]
    )

    plt.xlabel('# Predictions')
    plt.ylabel('Precision (%)')
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    #plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
    if out != None:
        plt.savefig(out)
    plt.legend().set_visible(False)
    plt.show()

def success_rate_IFSP(df, rank_t_list):
    results = []

    for origin_pred in df['origin_pred'].unique():
        subset = df[df['origin_pred'] == origin_pred] # subset of method predictions

        for rank_t in rank_t_list: # for each rank threshold
            rank_subset = subset[subset['RANK_pred'] <= rank_t]
            total_predictions = len(rank_subset)

            condition = rank_subset['TRUE'] == 1 # Correct predictions based on TRUE column
            correct_predictions = condition.sum()

            # Calculate the success rate
            success_rate = correct_predictions / total_predictions if total_predictions > 0 else 0

            # Calculate the standard error of the proportion
            se = (success_rate * (1 - success_rate) / total_predictions)**0.5 if total_predictions > 0 else 0

            # Calculate the confidence bounds
            lower_bound = max(0, success_rate - 1.96 * se)
            upper_bound = min(1, success_rate + 1.96 * se)

            results.append({
                'origin_pred': origin_pred,
                'rank_t': rank_t,
                'success_rate': round(success_rate, 3),
                'correct_predictions': correct_predictions,
                'total': total_predictions,
                'lower_se': round(lower_bound, 3),
                'upper_se': round(upper_bound, 3)
            })

    return pd.DataFrame(results)

# Outer loop to handle different distance thresholds
def calculate_success_rate_IFSP(df, dist_t_list, rank_t_list):
    combined_results = []

    for dist_t in dist_t_list:
        # Update the TRUE column based on the current distance threshold
        df['TRUE'] = (df['distance'] <= dist_t).astype(int)

        # Calculate success rate for the current distance threshold
        results_df = success_rate_IFSP(df, rank_t_list)
        results_df['dist_t'] = dist_t
        combined_results.append(results_df)

    return pd.concat(combined_results, ignore_index=True)

def success_rate_VNEGNN(df):
    results = []

    for origin_pred in df['origin_pred'].unique():
        subset = df[df['origin_pred'] == origin_pred] # subset of method predictions

        # Using RANK_pred <= N_SITES_ref
        rank_subset = subset[subset['RANK_pred'] <= subset['N_SITES_ref']]
        total_predictions = len(rank_subset)

        condition = rank_subset['TRUE'] == 1 # Correct predictions based on TRUE column
        correct_predictions = condition.sum()

        # Calculate the success rate
        success_rate = correct_predictions / total_predictions if total_predictions > 0 else 0

        # Calculate the standard error of the proportion
        se = (success_rate * (1 - success_rate) / total_predictions)**0.5 if total_predictions > 0 else 0

        # Calculate the confidence bounds
        lower_bound = max(0, success_rate - 1.96 * se)
        upper_bound = min(1, success_rate + 1.96 * se)

        results.append({
            'origin_pred': origin_pred,
            'success_rate': round(success_rate, 3),
            'correct_predictions': correct_predictions,
            'total': total_predictions,
            'lower_se': round(lower_bound, 3),
            'upper_se': round(upper_bound, 3)
        })

    return pd.DataFrame(results)

# Outer loop to handle different distance thresholds
def calculate_success_rate_VNEGNN(df, dist_t_list):
    combined_results = []

    for dist_t in dist_t_list:
        # Update the TRUE column based on the current distance threshold
        df['TRUE'] = (df['distance'] <= dist_t).astype(int)

        # Calculate success rate for the current distance threshold
        results_df = success_rate_VNEGNN(df)
        results_df['dist_t'] = dist_t
        combined_results.append(results_df)

    return pd.concat(combined_results, ignore_index=True)

def filter_redundant_pockets2(df, jaccard_threshold=0.75, distance_threshold=5.0):
    # Initialize an empty list to store the filtered rows
    filtered_rows = []

    # Group by 'origin' and 'rep_chain'
    try:
        grouped = df.groupby(['origin', 'rep_chain'])
    except:
        grouped = df.groupby(['origin_pred', 'rep_chain_ref'])

    for (origin, rep_chain), group in grouped:
        group = group.reset_index(drop=True)

        # Create a mask to keep track of rows to keep
        keep_mask = np.ones(len(group), dtype=bool)

        for i in range(len(group)):
            if not keep_mask[i]:
                continue
            try:
                up_aas_i = set(group.loc[i, 'up_aas'])
            except:
                up_aas_i = set(group.loc[i, 'up_aas_pred'])
            try:
                score_i = group.loc[i, 'SCORE']
            except:
                try:
                    score_i = group.loc[i, 'SCORE_pred']
                except:
                    score_i = group.loc[i, 'n_aas'] # for LIGYSIS, as there is no score, keeping larger site
            try:
                centre_i = np.array(group.loc[i, 'centre'])
            except:
                centre_i = np.array(group.loc[i, 'centre_pred'])

            for j in range(i + 1, len(group)):
                try:
                    up_aas_j = set(group.loc[j, 'up_aas'])
                except:
                    up_aas_j = set(group.loc[j, 'up_aas_pred'])
                try:
                    score_j = group.loc[j, 'SCORE']
                except:
                    try:
                        score_j = group.loc[j, 'SCORE_pred']
                    except:
                        score_j = group.loc[j, 'n_aas'] # for LIGYSIS, as there is no score, keeping larger site
                try:
                    centre_j = np.array(group.loc[j, 'centre'])
                except:
                    centre_j = np.array(group.loc[j, 'centre_pred'])

                jaccard = jaccard_index(up_aas_i, up_aas_j)
                distance = euclidean(centre_i, centre_j)

                if jaccard >= jaccard_threshold or distance <= distance_threshold:
                    if score_i >= score_j:
                        keep_mask[j] = False
                    else:
                        keep_mask[i] = False
                        break

        filtered_rows.extend(group[keep_mask].to_dict(orient='records'))

    return pd.DataFrame(filtered_rows)

def recalculate_IFSP_scores(df, base_dir):
    # Helper function to square and sum scores
    def square_and_sum(scores):
        return sum([x**2 for x in scores])

    # Initialize an empty DataFrame to collect results
    updated_dfs = []

    # Process each rep_chain separately
    for rep_chain in df['rep_chain'].unique():
        # Load the scores dictionary for the current rep_chain
        scores_dict_path = os.path.join(base_dir, rep_chain, f'{rep_chain}_lig_scores.pkl')
        scores_dict = read_from_pickle(scores_dict_path)

        # Filter the DataFrame for the current rep_chain
        rep_df = df[df['rep_chain'] == rep_chain]

        # Calculate the new SCORE by squaring the scores from the dictionary and summing them
        rep_df['SCORE'] = rep_df['aas'].apply(lambda aas_list: square_and_sum([scores_dict.get(aas, 0) for aas in aas_list]))

        # Sort by the new SCORE in descending order
        rep_df_sorted = rep_df.sort_values(by='SCORE', ascending=False)

        # Reset index to clean up index after sorting
        rep_df_sorted.reset_index(drop=True, inplace=True)

        # Reassign RANK based on new ordering
        rep_df_sorted['RANK'] = rep_df_sorted.index + 1

        # Collect the updated DataFrame
        updated_dfs.append(rep_df_sorted)

    # Concatenate all updated DataFrames from different rep_chains
    final_df = pd.concat(updated_dfs, ignore_index=True)

    return final_df

def success_rate_COMB(df, dist_t_list, rank_t_list, irel_t_list, TOTAL_preds = None):#, total_predictions=None):
    results = []

    for origin_pred in df['origin_pred'].unique():
        subset = df[df['origin_pred'] == origin_pred]  # subset of method predictions
        if TOTAL_preds != None:  # can pass absolute maximum (total observed pockets)
            total_predictions = TOTAL_preds
        else:
            total_predictions = len(subset)  # otherwise, this will be the number of pockets in rows for which the method predicts
        print(f"Total predictions for {origin_pred} is {total_predictions}")

        for dist_t in dist_t_list:  # for each distance threshold
            for rank_t in rank_t_list:  # for each rank threshold
                for irel_t in irel_t_list:  # for each intersection relative threshold
                    # Evaluate conditions with logical OR for distance and intersection relative thresholds
                    condition = ((subset['distance'] <= dist_t) | (subset['relative_intersection'] >= irel_t)) & \
                                (subset['RANK_pred'] <= (subset['N_SITES_ref'] + rank_t))
                    correct_predictions = condition.sum()

                    # Calculate the success rate
                    success_rate = correct_predictions / total_predictions if total_predictions > 0 else 0

                    # Calculate the standard error of the proportion
                    se = (success_rate * (1 - success_rate) / total_predictions)**0.5 if total_predictions > 0 else 0

                    # Calculate the confidence bounds
                    lower_bound = max(0, success_rate - 1.96 * se)
                    upper_bound = min(1, success_rate + 1.96 * se)

                    results.append({
                        'origin_pred': origin_pred,
                        'rank_t': rank_t,
                        'dist_t': dist_t,
                        'irel_t': irel_t,
                        'success_rate': round(success_rate, 3),
                        'correct_predictions': correct_predictions,
                        'total': total_predictions,
                        'lower_se': round(lower_bound, 3),
                        'upper_se': round(upper_bound, 3)
                    })

    return pd.DataFrame(results)

def plot_recall_heatmap(df):
    # Pivot the DataFrame to format suitable for heatmap
    heatmap_data = df.pivot_table(index='irel_t', columns='dist_t', values='success_rate', aggfunc=np.mean)
    
    # Plotting the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, cmap='inferno', annot=False, fmt=".2f", linewidths=.5, linecolor = "k")
    plt.title('Heatmap of Success Rate by irel_t and dist_t')
    plt.xlabel('dist_t')
    plt.ylabel('irel_t')
    plt.show()
    
def plot_results_def_NEW_2cols(df, methods, fix_col1, fix_col2, fix_val1, fix_val2, var_col, var_col_label, legend_order, ylabel="Recall", xline=None, xticks=None, xlim=None, yline=None, yticks=None, ylim=None, FSIZE=(10, 6), DPI=100, MS = 10, LW = 2.5, error = True, out = None):

    plt.figure(figsize=FSIZE, dpi=DPI)

    # Filter the DataFrame based on the two fixed column values
    df_filtered = df[(df[fix_col1] == fix_val1) & (df[fix_col2] == fix_val2)].query('origin_pred in @methods')

    df_filtered["success_rate"] = df_filtered["success_rate"]*100
    df_filtered["lower_se"] = df_filtered["lower_se"]*100
    df_filtered["upper_se"] = df_filtered["upper_se"]*100

    # Determine linestyle based on the presence of 'PRANK' in the 'origin_pred' column
    #linestyles = {origin_pred: (1, 1) if 'PRANK' in origin_pred else '' for origin_pred in df_filtered['origin_pred'].unique()}

    if var_col == "rank_t":
        df_filtered.loc[df['rank_t'] > 20, 'rank_t'] = 22
    
    sns.lineplot(data=df_filtered, x=var_col, y='success_rate', hue='origin_pred',
                 style='origin_pred', markers=markers_dict, palette=palette, 
                 markersize=MS, linewidth=LW, err_style=None, hue_order=legend_order, style_order=legend_order, 
                 dashes=linestyles_dict)

    if error:
        for name, group in df_filtered.groupby('origin_pred'):
            plt.errorbar(group[var_col], group['success_rate'],
                         yerr=[group['success_rate'] - group['lower_se'], group['upper_se'] - group['success_rate']],
                         fmt='none', ecolor=palette[name], capsize=3)

    plt.xlabel(var_col_label)
    plt.ylabel(ylabel)

    plt.legend().set_visible(False)
    
        

    if xticks is not None:
        if var_col != 'rank_t':
            plt.xticks(xticks)
        else:
            new_xticks = list(xticks)+[22]
            new_labs = [f'n+{str(el)}' for el in xticks] + ["ALL"]
            plt.xticks(new_xticks, new_labs, rotation = -45)

    if xlim is not None:
        plt.xlim(xlim)

    if xline is not None:
        plt.axvline(x=xline, linestyle="--", color="k", linewidth=1)

    if yticks is not None:
        plt.yticks(yticks)

    if ylim is not None:
        plt.ylim(ylim)

    if yline is not None:
        plt.axhline(y=yline, linestyle="--", color="k", linewidth=1)

    if out != None:
        plt.savefig(out)

    plt.show()

def normalise_scores(df, methods):
    # Work on a copy of the DataFrame to avoid modifying the original data
    result_df = df.copy()
    result_df['score'] = result_df['score'].clip(upper=1)
    return result_df
    
### VARIABLES ###

palette = {
    'LIGYSIS': 'forestgreen',
    'VN-EGNN': 'orchid',
    'VN-EGNN-NR': 'orchid',
    'IF-SitePred': 'dodgerblue',
    'IF-SitePred-NR': 'dodgerblue',
    'IF-SitePred-rescored-NR': 'dodgerblue',
    'GrASP': 'firebrick',
    'PUResNet': 'royalblue',
    'PUResNet-AA': 'royalblue',
    'PUResNet+PRANK': 'royalblue',
    'PUResNet+PRANK+Cons': 'royalblue',
    'DeepPocket-Segmented': 'sienna',
    'DeepPocket-Segmented-NR': 'sienna',
    'DeepPocket-Rescored': 'tan',
    'P2Rank+Cons': 'magenta',
    'P2Rank': 'orange',
    'fpocket+PRANK+Cons': 'gray',
    'fpocket+PRANK': 'black',
    'fpocket': 'gray',
    'PocketFinder': 'blue',
    'PocketFinder-AA': 'blue',
    'PocketFinder-SS': 'blue',
    'PocketFinder+PRANK': 'blue',
    'PocketFinder+PRANK+Cons': 'blue',
    'Ligsite': 'red',
    'Ligsite-AA':  'red',
    'Ligsite-SS': 'red',
    'Ligsite+PRANK': 'red',
    'Ligsite+PRANK+Cons': 'red',
    'Surfnet': 'green',
    'Surfnet-AA': 'green',
    'Surfnet-SS': 'green',
    'Surfnet+PRANK': 'green',
    'Surfnet+PRANK+Cons': 'green',
 }

markers_dict = {
    "VN-EGNN": '8',
    "VN-EGNN-NR": '8',
    "IF-SitePred": 'D',
    "IF-SitePred-NR": 'D',
    "IF-SitePred-rescored-NR": 'D',
    "GrASP": 's',
    "PUResNet": 'p',
    "PUResNet-AA": 'p', 
    "PUResNet+PRANK": 'p',
    "PUResNet+PRANK+Cons": 'p',
    "DeepPocket-Segmented": 'H',
    "DeepPocket-Segmented-NR": 'H',
    "DeepPocket-Rescored": 'h',
    "P2Rank+Cons": 'v',
    "P2Rank": '^',
    'fpocket+PRANK+Cons': '.',
    "fpocket+PRANK": '.',
    "fpocket": 'o',
    "PocketFinder": '*',
    "PocketFinder-AA": "*",
    "PocketFinder-SS": '*',
    "PocketFinder+PRANK": '*',
    "PocketFinder+PRANK+Cons": '*',
    "Ligsite": 'P',
    "Ligsite-AA": 'P',
    "Ligsite-SS": 'P',
    "Ligsite+PRANK": 'P',
    "Ligsite+PRANK+Cons": 'P',
    "Surfnet": 'X',
    "Surfnet-AA": 'X',
    'Surfnet-SS': 'X',
    "Surfnet+PRANK": 'X',
    "Surfnet+PRANK+Cons": 'X',
}

linestyles_dict = {
    "VN-EGNN": '',
    "VN-EGNN-NR": (1, 1),                    ## keeping this one
 
    "IF-SitePred": '',
    "IF-SitePred-NR": (1, 1),
    "IF-SitePred-rescored-NR": (5, 1),       ## keeping this one

    "GrASP": '',                             ## no variants

    "PUResNet": '',
    "PUResNet-AA": (3, 1, 1, 1, 1, 1),
    "PUResNet+PRANK": (3, 1, 1, 1),          ## keeping this one
    "PUResNet+PRANK+Cons": (1, 1),

    "DeepPocket-Segmented": '', 
    "DeepPocket-Segmented-NR": (1, 1),       ## keeping this one

    "DeepPocket-Rescored": '',               ## no variants
    "P2Rank+Cons": '',                       ## no variants
    "P2Rank": '',                            ## no variants
    "fpocket+PRANK+Cons": (1, 1),
    #"fpocket+PRANK": (3, 1, 1, 1),           ## no variants
    "fpocket+PRANK": '',                     ## no variants
    "fpocket": '',                           ## no variants

    "PocketFinder": '',
    "PocketFinder-AA": (3, 1, 1, 1, 1, 1),
    "PocketFinder-SS": (5, 10),              ## keeping this one
    "PocketFinder+PRANK": (3, 1, 1, 1),
    "PocketFinder+PRANK+Cons": (1, 1),

    "Ligsite": '',
    "Ligsite-AA": (3, 1, 1, 1, 1, 1),
    "Ligsite-SS": (5, 10),                   ## keeping this one
    "Ligsite+PRANK": (3, 1, 1, 1),
    "Ligsite+PRANK+Cons": (1, 1),

    "Surfnet": '',
    "Surfnet-AA": (3, 1, 1, 1, 1, 1),
    'Surfnet-SS': (5, 10),                   ## keeping this one
    "Surfnet+PRANK": (3, 1, 1, 1),
    "Surfnet+PRANK+Cons": (1, 1),
}

linestyles_dict2 = {
    "VN-EGNN": (0, ()),
    "VN-EGNN-NR": (0, (2, 1)),                    ## keeping this one
 
    "IF-SitePred": (0, ()),
    "IF-SitePred-NR": (0, (2, 1)),
    "IF-SitePred-rescored-NR": (0, (5, 1)),       ## keeping this one

    "GrASP": (0, ()),                             ## no variants

    "PUResNet": (0, ()),
    "PUResNet-AA": (0, (5, 1)),
    "PUResNet+PRANK": (0, (3, 1, 1, 1)),          ## keeping this one
    "PUResNet+PRANK+Cons": (0, (1, 1)),

    "DeepPocket-Segmented": (0, ()),
    "DeepPocket-Segmented-NR":(0, (2, 1)),       ## keeping this one

    "DeepPocket-Rescored": (0, ()),               ## no variants
    "P2Rank+Cons": (0, ()),                       ## no variants
    "P2Rank": (0, ()),                            ## no variants

    "fpocket+PRANK+Cons": (0, (1, 1)),
    "fpocket+PRANK": (0, (3, 1, 1, 1)),           ## no variants
    "fpocket+PRANK": (0, ()),                     ## no variants
    "fpocket": (0, ()),                           ## no variants

    "PocketFinder": (0, ()),
    "PocketFinder-AA": (0, (5, 1)),
    'PocketFinder-SS': (0, (5, 10)),                   ## keeping this one
    "PocketFinder+PRANK": (0, (3, 1, 1, 1)),
    "PocketFinder+PRANK+Cons": (0, (1, 1)), 

    "Ligsite": (0, ()),
    "Ligsite-AA": (0, (5, 1)),
    'Ligsite-SS': (0, (5, 10)),                   ## keeping this one
    "Ligsite+PRANK": (0, (3, 1, 1, 1)),
    "Ligsite+PRANK+Cons": (0, (1, 1)), 

    "Surfnet": (0, ()),
    "Surfnet-AA": (0, (5, 1)),
    'Surfnet-SS': (0, (5, 10)),                   ## keeping this one
    "Surfnet+PRANK": (0, (3, 1, 1, 1)),
    "Surfnet+PRANK+Cons": (0, (1, 1)), 
}


volumes_dict= read_from_pickle("./results/MASTER_POCKET_SHAPE_DICT_EXTENDED_TRANS.pkl")
