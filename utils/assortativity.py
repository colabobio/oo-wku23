import numpy as np
from statsmodels.stats.weightstats import DescrStatsW

def get_attrib_data(G, attrib_name):
    data = []
    for n, m, d in G.edges(data=True):
        rn = G.nodes[n][attrib_name]
        rm = G.nodes[m][attrib_name]
        w = d['weight']
        if not np.isnan(rn) and not np.isnan(rm):
            data.append([rn, rm, w])
    return np.array(data)

def calc_weighted_attrib_assortativity(G, attrib):
    data = get_attrib_data(G, attrib)
    wdata = DescrStatsW(data[:,0:2], weights=data[:,2])
    return wdata.corrcoef[0, 1]

def calc_weighted_attrib_assortativity_sig(G, attrib, rho, n_boot=1000):
    data = get_attrib_data(G, attrib)
    p_value = 0.0
    for bootstrap in range(n_boot):
        rho_boot = 0
        xy = np.random.choice(data[:,0:2].flatten(), size=(data.shape[0], 2))
        ww = np.random.choice(data[:,2], size=data.shape[0])
        wdata = DescrStatsW(xy, weights=ww)
        rho_boot = wdata.corrcoef[0, 1]
        if rho_boot >= rho:
            p_value += 1
    p_value /= n_boot
    return p_value

def calc_weighted_assortativity(data, calc_pval=False, n_boot=100):
    W = np.sum(data[:,2])
    
    Xbar = 0
    Ybar = 0
    Xsig = 0
    Ysig = 0
    for row in range(data.shape[0]):
        Xbar += data[row,0]*data[row,2]/W
        Ybar += data[row,1]*data[row,2]/W
    for row in range(data.shape[0]):
        Xsig += data[row,2]*(data[row,0]-Xbar)**2/W
        Ysig += data[row,2]*(data[row,1]-Ybar)**2/W
    Xsig = np.sqrt(Xsig)
    Ysig = np.sqrt(Ysig)

    rho = 0
    for row in range(data.shape[0]):
        rho +=data[row,2]*(data[row,0]-Xbar)*(data[row,1]-Ybar)
        rho /= W*Xsig*Ysig

    if calc_pval:
        p_value = 0.0
        Xbar_boot = np.mean(data[:,0])
        Xsig_boot = np.std(data[:,0])
        for bootstrap in range(n_boot):
            rho_boot = 0
            for row in range(data.shape[0]):
                x = np.random.choice(data[:,0])
                y = np.random.choice(data[:,1])
                w = np.random.choice(data[:,2])
                rho_boot += w*(x-Xbar_boot)*(y-Xbar_boot)
                rho_boot /= W*Xsig_boot**2
            if rho_boot >= rho:
                p_value += 1
                # print(p_value/(bootstrap+1))
        p_value /= n_boot
        return rho, p_value
    
    return rho, None