import pandas as pd
import numpy as np
from scipy.stats import poisson
import statsmodels.api as sm
import networkx as nx

# Returns super-spreader for a given tree
def ss_stats_per_tree(tree, ss_threshold, root, id):
    # Convert tree node labels to string, assuming they might not be
    tree = nx.relabel_nodes(tree, {node: str(node) for node in tree.nodes()})
    root = str(root)  # Ensure root is a string to match node labels

    names_ss = [node for node, out_degree in tree.out_degree() if out_degree > ss_threshold]

    ss_success_rate_num = ss_success_rate_denom = case_success_rate_num = case_success_rate_denom = 0
    infectors_ris = []
    non_infectors_ris = []

    # first check if tree has superspreaders
    if len(names_ss) == 0:
         # success rates irrelevant here since no superspreaders
        case_success_rate_num = case_success_rate_denom = ss_success_rate_num = ss_success_rate_denom = None
        for node in tree.nodes:
            non_infectors_ris.append(tree.out_degree(node))
    else:
        # if the tree does have superspreaders:
        names_dir_infected = [] # list of nodes directly infected by ALL superspreaders
        infectors = []          # list of nodes that infected superspreaders 

        for ss in names_ss: # for each superspreader
            # add names of those directly infected by ss to list
            directly_infected = [n for n in tree.neighbors(ss)]
            names_dir_infected.extend(directly_infected)

            # Find the infector of the superspreader, if any
            infector = [node for node in tree.predecessors(ss)]
            if infector:
                # add all infectors of superspreaders in this tree to this list
                infectors.extend(infector)

        # all ppl infected by ss, moved to not inside loop
        ss_success_rate_denom = len(set(names_dir_infected))
        ss_success_rate_num = len([name for name in names_ss if name in names_dir_infected]) # if ss infected by another ss

        # for each node, see if infected by a ss, calculate its out degree, also calculate success rate nums and denoms
        for node in tree.nodes:
            # if is a superspreader
            if node in names_ss:
                if node in infectors: # if infected another superspreader
                    infectors_ris.append(tree.out_degree(node))
                else: # if did not infect another superspreader
                    non_infectors_ris.append(tree.out_degree(node))
            else:
                case_success_rate_denom += tree.out_degree(node)
                if node in infectors: # if infected superspreader
                    case_success_rate_num += 1
                    infectors_ris.append(tree.out_degree(node))
                else: # if did not infect superspreader
                    non_infectors_ris.append(tree.out_degree(node))

    # Calculate success rates, taking care of division by zero
    ss_success_rate = ss_success_rate_num / ss_success_rate_denom if ss_success_rate_denom else None
    case_success_rate = case_success_rate_num / case_success_rate_denom if case_success_rate_denom else None

    # Calculate superspreader dyads metrics
    num_ss = len(names_ss)
    size = tree.number_of_nodes()
    num_terminal = sum(1 for _, out_degree in tree.out_degree() if out_degree == 0)
    observed_ss_dyads = ss_success_rate_num if ss_success_rate_num is not None else 0
    expected_ss_dyads = (num_ss * (num_ss - 1)) / size
    expected_ss_dyads_nt = (num_ss * (num_ss - 1)) / (size - num_terminal) if (size - num_terminal) else 0
    excess_ss_dyads = observed_ss_dyads - expected_ss_dyads
    excess_ss_dyads_nt = observed_ss_dyads - expected_ss_dyads_nt

    output = {
        "SS Success Rate": ss_success_rate,
        "Case Success Rate": case_success_rate,
        "Obs SS Dyads": observed_ss_dyads,
        "Exp SS Dyads": expected_ss_dyads,
        "Excess SS Dyads": excess_ss_dyads,
        "Exp SS Dyads JLS": expected_ss_dyads_nt,
        "Excess SS Dyads JLS": excess_ss_dyads_nt,
        "Infectors Ris": infectors_ris,
        "Non-Infectors Ris": non_infectors_ris
    }

    return output

# Determines average out degree of nodes infected vs not infected by a superspreader, as well as success rate
def inf_by_ss_info(tree, thresh, thresh_nt, thresh_nlt, thresh_g01, thresh_g01_nt):
    # Ensure the tree is correctly prepared for processing
    #tree = nx.relabel_nodes(tree, {node: str(node) for node in tree.nodes()})
    root = [n for n, d in tree.in_degree() if d == 0][0]  # Identifying the root node
    
    # Main statistics calculation for different thresholds
    all_nodes = ss_stats_per_tree(tree, thresh, root, id)
    noterm_nodes = ss_stats_per_tree(tree, thresh_nt, root, id)
    nolastterm_nodes = ss_stats_per_tree(tree, thresh_nlt, root, id)
    g01_nodes = ss_stats_per_tree(tree, thresh_g01, root, id)
    noterm_g01_nodes = ss_stats_per_tree(tree, thresh_g01_nt, root, id)

    # Packaging the results into a single dictionary
    bound_output = {
        "all_nodes": all_nodes,
        "noterm_nodes": noterm_nodes,
        "nolastterm_nodes": nolastterm_nodes,
        "g01_nodes": g01_nodes,
        "noterm_g01_nodes": noterm_g01_nodes
    }

    return bound_output

# Gives the ss stats for the first and second halves of the tree
def gen_halves_stats(tree, thresh, n_reps=10):
    root = [n for n, d in tree.in_degree() if d == 0][0]
    gens = np.max(list(nx.single_source_shortest_path_length(tree, root).values()))
    output_all = []

    # do this nrepeats times to account for stochasticity
    for rep in range(n_reps):
        distances = nx.single_source_shortest_path_length(tree, root)
        first_nodes = [node for node, distance in distances.items() if distance < gens / 2]
        second_nodes = [node for node, distance in distances.items() if distance > gens / 2]
        middle_nodes = [node for node, distance in distances.items() if distance == gens / 2]

        # randomly assign nodes to first and second halves
        if middle_nodes:
            random_nums = np.random.uniform(size=len(middle_nodes))
            for i, node in enumerate(middle_nodes):
                if random_nums[i] <= 0.5:
                    first_nodes.append(node)
                else:
                    second_nodes.append(node)

        # Calculate statistics for first and second halves
        def calc_stats(nodes):            
            out_degrees = [tree.out_degree(node) for node in nodes]
            out_degrees_nt = [deg for deg in out_degrees if deg > 0]

            # all
            r_a = np.mean(out_degrees)
            thresh_a = poisson.ppf(0.99, r_a)
            dp_a = get_nbinom_size(fit_nbinom(out_degrees))
            num_ss_a = len([deg for deg in out_degrees if deg > thresh])
            ss_freq_a = num_ss_a / len(nodes) if nodes else 0
            
            # no term
            r_nt = np.mean(out_degrees_nt)
            thresh_nt = poisson.ppf(0.99, r_nt)
            dp_nt = get_nbinom_size(fit_nbinom(out_degrees_nt))
            num_ss_nt = len([deg for deg in out_degrees_nt if deg > thresh_nt])
            ss_freq_nt = num_ss_nt / len(nodes) if nodes else 0
            
            return r_a, dp_a, ss_freq_a, r_nt, dp_nt, ss_freq_nt

        r_1_a, dp_1_a, ss_freq_1_a, r_1_nt, dp_1_nt, ss_freq_1_nt = calc_stats(first_nodes)
        r_2_a, dp_2_a, ss_freq_2_a, r_2_nt, dp_2_nt, ss_freq_2_nt = calc_stats(second_nodes)

        output = pd.DataFrame({
            'Gens': [gens],
            'R_1': [r_1_a],
            'R_2': [r_2_a],
            'DP_1': [dp_1_a],
            'DP_2': [dp_2_a],
            'SS Freq 1': [ss_freq_1_a],
            'SS Freq 2': [ss_freq_2_a],            
            'R_1_NT': [r_1_nt],
            'R_2_NT': [r_2_nt],
            'DP_1_NT': [dp_1_nt],
            'DP_2_NT': [dp_2_nt],
            'SS Freq 1 NT': [ss_freq_1_nt],
            'SS Freq 2 NT': [ss_freq_2_nt],
            'rep_num': [rep]
        })

        output_all.append(output)

    output_all_df = pd.concat(output_all, ignore_index=True)
    summary_stats = output_all_df.mean().to_frame().T

    return {"summary": summary_stats, "indiv_sims": output_all_df}

def extract_ss_info_all(tree, ss_threshold):
    size = tree.number_of_nodes()
    root = [n for n, d in tree.in_degree() if d == 0][0]  # Identifying the root node
    names_ss = [n for n, d in tree.out_degree() if d > ss_threshold]  # Superspreaders

    # Check if there are superspreaders
    if names_ss:
        num_dir_infected = []
        num_down_infected = []
        ss_dir_infect = []
        num_infector_infected = []

        for ss in names_ss:
            # Directly infected
            num_dir_infected.append(tree.out_degree(ss))
            # Indirectly infected (downstream)
            num_down_infected.append(len(nx.descendants(tree, ss)))
            
            # Infector analysis
            predecessors = list(tree.predecessors(ss))
            if predecessors:
                infector = predecessors[0]
                num_infector_infected.append(tree.out_degree(infector))
                ss_dir_infect.append(infector in names_ss)
            else:
                num_infector_infected.append(None)
                ss_dir_infect.append(None)

        max_down = np.max(num_down_infected)
        output = pd.DataFrame({
            'SS': names_ss,
            'Thresh': ss_threshold,
            'Direct Infected': num_dir_infected,
            'Direct Prop': [n / (size - 1) for n in num_dir_infected],
            'Indiv Downstream': num_down_infected,
            'Indiv Downstream Prop': [n / (size - 1) for n in num_down_infected],
            'Max Downstream Prop': [max_down / (size - 1)] * len(names_ss),
            'Num SS in Tree': [len(names_ss)] * len(names_ss),
            'Directly Infected by SS': ss_dir_infect,
            'Size': [size] * len(names_ss),
            'R Infector': num_infector_infected,
        })
    else:
        output = pd.DataFrame()

    return output

# Utility functions

def fit_nbinom0(data):
    """Fit a negative binomial distribution to data and return the size parameter."""
    mean = np.mean(data)
    variance = np.var(data)
    if variance > mean:
        # Estimate the parameters for the negative binomial
        p = mean / variance
        r = mean * p / (1 - p)
        return {"size": r}
    else:
        return {"size": np.nan}

def fit_nbinom(data, summary=False):
    # https://anton-granik.medium.com/fitting-and-visualizing-a-negative-binomial-distribution-in-python-3cc27fbc7ecf
    res = sm.NegativeBinomial(data, np.ones_like(data)).fit(start_params=[1,1])
    if summary: print(res.summary())
    return res.params

def get_nbinom_size(params):
    mu = np.exp(params[0])
    p = 1 / (1 + mu * params[1])
    n = np.min([100, mu * p/(1-p)]) # NOTE: The fitting does not work for second half degrees... 
    return n

def plot_nbinom_fit(data, params):
    mu = np.exp(params[0])
    p = 1 / (1 + mu * params[1])
    n = mu * p/(1-p)
    x = np.linspace(0, max(data), max(data)+1)    
    ax = sns.distplot(data, kde = False, norm_hist=True, label='Real values')
    ax.plot(x, nbinom.pmf(x, n, p), 'g-', lw=2, label='Fitted NB')
    leg = ax.legend()
    plt.title('Real vs Fitted NB Dist')