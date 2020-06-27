import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def extract(file):
    if file.endswith('.vtp'):
        import base64
        from xml.dom.minidom import parse
        def to_np(base64_array):
            return np.frombuffer(base64.decodebytes(base64_array.encode("ascii")), dtype=np.dtype(np.float32).newbyteorder("<"))
        def get_data(document):
            arrays = document \
            .getElementsByTagName("VTKFile")[0] \
            .getElementsByTagName("PolyData")[0] \
            .getElementsByTagName("Piece")[0] \
            .getElementsByTagName("PointData")[0] \
            .getElementsByTagName("DataArray")
            solution = next(filter(lambda a: a.getAttribute("Name") == "solution", arrays))
            return solution.childNodes[0].data
        return to_np(get_data(parse(file)))[1:]
    if file.endswith('.py'):
        import sys
        if 'py_reader' not in sys.modules:
            sys.path.append("/home/huberfx/opendihu/scripts")
            import py_reader
        data = py_reader.load_data([file])
        data = data[0]['data']
        solution = next(filter(lambda d: d['name'] == 'solution', data))
        component0 = next(filter(lambda d: d['name'] == '0', solution['components']))
        return component0['values']
    raise "FileType not understood: "+file
        
def extract4(file):
    if file.endswith('.vtp'):
        a = extract(file)
        # hack: last datapoint is broken?
        # a = a[:a.shape[0]//4*4]
        # a = a[1:]
        return a.reshape((-1,4)) # a[x, component]
    if file.endswith('.py'):
        import sys
        if 'py_reader' not in sys.modules:
            sys.path.append("/home/huberfx/opendihu/scripts")
            import py_reader
        data = py_reader.load_data([file])
        data = data[0]['data']
        solution  = next(filter(lambda d: d['name'] == 'solution', data))
        componentX = lambda x: next(filter(lambda d: d['name'] == str(x), solution['components']))
        return np.vstack([componentX(i)['values'] for i in "Vmhn"]).T
    raise "FileType not understood: "+file






def opt_coarsen(data, opt_start=(0,0), opt_0 = (0,0), eps = np.infty, norm=2):
    def coarsen_N(ix):
        return (ix[0] + 1, ix[1])
    def coarsen_M(ix):
        return (ix[0], ix[1] + 1)
    dat_0 = data[opt_0]
    opt_curr = opt_start
    curr  = data[opt_curr]
    opts = []
    opts.append(opt_curr)
    try:
        while np.linalg.norm(curr - dat_0, ord=norm) < eps:
            c_N  = data[coarsen_N(opt_curr)]
            c_M  = data[coarsen_M(opt_curr)]
            e_N  = np.linalg.norm(c_N - dat_0, ord=norm)
            e_M  = np.linalg.norm(c_M - dat_0, ord=norm)
            if e_N > e_M:
                opt_curr = coarsen_M(opt_curr)
            else:
                opt_curr = coarsen_N(opt_curr)
            curr = data[opt_curr]
            opts.append(opt_curr)
    except IndexError as e:
        print(e)
    return np.array(opts)

def opt_refine(data, opt_start=(-1,-1), eps = np.infty, norm=2):
    if opt_start[0] < 0: opt_start = (data.shape[0] + opt_start[0], opt_start[1])
    if opt_start[1] < 0: opt_start = (opt_start[0], data.shape[1] + opt_start[1])
    def refine_N(ix):
        return (ix[0] - 1, ix[1])
    def refine_M(ix):
        return (ix[0], ix[1] - 1)
    opt_curr = opt_start
    curr  = data[opt_curr]
    opts = []
    c_NM = curr
    opts.append(opt_curr)
    try:
        while np.linalg.norm(curr - c_NM, ord=norm) < eps and opt_curr[0] > 0 and opt_curr[1] > 0:
            c_N  = data[refine_N(opt_curr)]
            c_M  = data[refine_M(opt_curr)]
            c_NM  = data[refine_M(refine_N(opt_curr))]
            e_N  = np.linalg.norm(c_N - c_NM, ord=norm)
            e_M  = np.linalg.norm(c_M - c_NM, ord=norm)
            e_NM  = np.linalg.norm(c_NM - c_NM, ord=norm)
            if e_N > e_M:
                opt_curr = refine_M(opt_curr)
            else:
                opt_curr = refine_N(opt_curr)
            curr = data[opt_curr]
            opts.append(opt_curr)
    except IndexError as e:
        print(e)
    return np.array(opts)

# time = 42
# MS = [1,2,3,4,5,6,7,8,9,10,16,32,64]
# MS = [2,4,6,8,10,16,32,64,128,256]
# data = np.array([extract("out/fibre_dt_{}_{:07}.vtp".format(M,time)) for M in MS])
# data = np.array([extract("out_neon_IE/fibre_dt_{}_{}_{}_{:07}.vtp".format("3e-3",1,M,time)) for M in MS])
# surplus = data[1:] - data[:-1]
# err = data[:] - data[-1]
# norms_sur_2 = np.linalg.norm(surplus,axis=1)
# norms_sur_oo = np.linalg.norm(surplus,ord=np.inf,axis=1)
# norms_err_2 = np.linalg.norm(err,axis=1)
# norms_err_oo = np.linalg.norm(err,ord=np.inf,axis=1)
#
# plt.subplot(231)
# plt.plot(surplus.T)
# plt.legend(MS)
# plt.subplot(232)
# plt.semilogy(np.abs(surplus.T))
# plt.subplot(234)
# plt.plot(err.T)
# plt.subplot(235)
# plt.semilogy(np.abs(err.T))
# plt.subplot(233)
# plt.loglog(MS[:-1], norms_sur_2, "-o", basey=2, basex=2)
# plt.loglog(MS[:-1], norms_sur_oo, "-o", basey=2, basex=2)
# plt.legend(["2", "oo"])
# plt.subplot(236)
# plt.loglog(MS[:-1], norms_err_2[:-1], "-o", basey=2, basex=2)
# plt.loglog(MS[:-1], norms_err_oo[:-1], "-o", basey=2, basex=2)
# plt.show()
#
# print(norms_sur_2)
# print(norms_sur_oo)
# print(norms_err_2)
# print(norms_err_oo)

def nanan(x):
    x = np.array(x)
    # return x
    x[np.isinf(x)] = np.min(x[np.logical_not(np.isinf(x))]) - 1
    return x

def ID(exp):
    return {1:"6e-3", 0:"3e-3"}.get(exp, "3e-3*2**{}".format(exp))


time = 1
exp = -1
DTS="6e-3"
DTS=ID(exp)
MS=[256,128,64,32,16,10,8,6,4,2] # dt_3D/M for reaction
MS=[256,128,64,32,16,8,4,2] # dt_3D/M for reaction
NS=[32,16,8,4,3,2,1] # dt_3D/N for diffusion
NS=[32,16,8,4,2,1] # dt_3D/N for diffusion

MS=[256,128,64,32,16,8,4,2] # dt_3D/M for reaction
NS=[256,128,64,32,16,8,4,2,1] # dt_3D/N for diffusion

exp_sol = -3
time_sol = time*2**(exp-exp_sol)
print("Exp:", exp, "exp sol:", exp_sol)
print("Time:", time, "time sol:", time_sol)
sol_r2 = extract4("out_neon_CN_magic_480/fibre_dt_r_short_{}_{}_{}_{:07}.vtp".format(ID(exp_sol),1,2,time_sol))
data_d = np.array([[extract("out_neon_CN_magic_480/fibre_dt_short_{}_{}_{}_{:07}.vtp".format(DTS,N,M,time)) for N in NS] for M in MS])
data_r2 = np.array([[extract4("out_neon_CN_magic_480/fibre_dt_r_short_{}_{}_{}_{:07}.vtp".format(DTS,N,M,time)) for N in NS] for M in MS])
print(sol_r2.shape)
print(data_d.shape)
print(data_r2.shape)
print("was ist mit den shapes los? (in dem python callback sind es 1191 nodes)")


# for ix in range(1):
for ix in range(1):
    data_s = data_d[:,:,:]
    sol_s = sol_r2[:]
    data_s = data_r2[:,:,:,ix]
    sol_s = sol_r2[:,ix]
    plt.figure("Err After Splitting and Diffusion. Component {}".format(ix))

    def err(data):
        return data - data[0,0,:]
    def err_true(data):
        return data - sol_s
    def norm_2(data):
        return np.linalg.norm(data, axis=2)
    def norm_oo(data):
        return np.linalg.norm(data, ord=np.inf, axis=2)
    def surp_2(data):
        return np.linalg.norm(data[1:,1:] - data[1:,:-1] - data[:-1,1:] + data[:-1,:-1], axis=2)
    def minmax(data):
        return np.min(data), np.max(data)

    err_s = err(data_s)
    err_s_2 = norm_2(err_s)
    err_s_oo = norm_oo(err_s)
    surp_s_2 = surp_2(data_s)
    err_d = err_true(data_s)
    err_d_2 = norm_2(err_d)
    err_d_oo = norm_oo(err_d)
    # surp_d_2 = surp_2(sol_s)

    mi,ma = minmax(np.array([minmax(nanan(np.log(np.abs(err_s_2)))), minmax(nanan(np.log(np.abs(err_d_2))))]))
    mima2 = {'vmin':mi, 'vmax':ma}
    print(mima2)

    plt.subplot(231)
    plt.imshow(nanan(np.log(np.abs(err_s_2))), cmap="gray_r", origin="lower", **mima2)
    plt.title("||.||_2")
    plt.xlabel("N (diffusion)")
    plt.ylabel("M (reaction)")
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_s, (0,0)).T[::-1])
    plt.plot(*opt_coarsen(data_s, (1,1)).T[::-1] + 0.1)
    plt.subplot(232)
    plt.title("||.||_oo")
    plt.imshow(nanan(np.log(np.abs(err_s_oo))), cmap="gray_r", origin="lower")
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_s, (0,0), norm=np.inf).T[::-1])
    plt.plot(*opt_coarsen(data_s, (1,1), norm=np.inf).T[::-1] + 0.1)
    plt.subplot(233)
    plt.title("surplus")
    try:
        plt.imshow(nanan(np.log(np.abs(surp_s_2))), cmap="gray", origin="lower")
    except: pass
    plt.xticks(np.arange(len(NS)-1), NS)
    plt.yticks(np.arange(len(MS)-1), MS)
    plt.colorbar()

    plt.subplot(234)
    plt.imshow(nanan(np.log(np.abs(err_d_2))), cmap="gray_r", origin="lower", **mima2)
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.xlabel("N (diffusion)")
    plt.ylabel("M (reaction)")
    plt.colorbar()
    plt.subplot(235)
    plt.imshow(nanan(np.log(np.abs(err_d_oo))), cmap="gray_r", origin="lower")
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.colorbar()
    # plt.subplot(236)
    # plt.imshow(nanan(np.log(np.abs(surp_d_2))), cmap="gray", origin="lower")
    # plt.xticks(np.arange(len(NS)-1), NS)
    # plt.yticks(np.arange(len(MS)-1), MS)
    # plt.colorbar()

    # plt.show(block=False)

    # mi,ma = minmax(err(data_s))
    # mi,ma = (None,None)
    plt.figure("Errors After Splitting. Component {}. rows = M (reaction) cols = N (diffusion).".format(ix))

    max_e = np.max(np.abs(err_d))
    max_sol = np.max(np.abs(data_s[0,0]))
    for ni in range(len(NS)):
        for mi in range(len(MS)):
            # nrows ncols index (topleft to right)
            plt.subplot(len(MS), len(NS), 1 + ni + mi*len(NS))
            # plt.ylim(mi,ma)
            plt.axvline(np.argmax(np.abs(data_s[len(MS) - 1 - mi,ni,:data_s.shape[2]//2])), color="red")
            plt.axvline(data_s.shape[2]//2+np.argmax(np.abs(data_s[len(MS) - 1 - mi,ni,data_s.shape[2]//2:])), color="red")
            # plt.axvline(np.argmax(np.abs(err(data_s)[len(MS) - 1 - mi,ni])), color="green", alpha=1)
            e = err_d[len(MS) - 1 - mi,ni]
            plt.xticks([])
            plt.yticks([])
            if False:
                plt.semilogy(np.abs(np.fft.fft(e))) # high freq in center low freq at borders
            else:
                plt.ylim(None, max_e)
                # plt.semilogy(np.abs(e))
                plt.plot(e)
                plt.twinx().plot(data_s[0,0], color="black", alpha=0.25)
                
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
                plt.twinx().tick_params(axis='y', which='both', left=False, right=False, labelright=False)
    # plt.show(block=False)











factor_NM = (1, 1)
factor_NM = (4, 1)
EXPS = [-3, -2, -1, 0, 1]
NMS = None
NS=[256,128,64,32,16,8,4,2,1,2**-1,2**-2,2**-3,2**-4]
MS=[256,128,64,32,16,8,4,2,1,2**-1,2**-2,2**-3,2**-4]
if factor_NM[0] == 1:
     NMS = [(N, factor_NM[1]*N) for N in NS if factor_NM[1]*N in MS]
elif factor_NM[1] == 1:
    NMS = [(factor_NM[0]*M, M) for M in MS if factor_NM[0]*M in NS]

# factor so that all methods stop after one time step on the coarsest level
def time_factor(exp, base=EXPS[-1]):
    return 2**(base-exp)

# factor to scale N and m so that all methods use the same time step
def NM_factor(exp, base=EXPS[-1]):
    # return 1
    return time_factor(EXPS[0],EXPS[-1])//time_factor(exp,base)

def try_load(file):
    try:
        return extract4(file)
    except FileNotFoundError:
        print("No file", file)
        return np.zeros(shape=(1191,4))*np.nan
# print ([["out_neon_CN_magic_480/fibre_dt_short_{}_{}_{}_{:07}.vtp".format(ID(E),N,M,time*time_factor(E)) for N,M in NMS] for E in EXPS])
data_e = np.array([[try_load("out_neon_CN_magic_480/fibre_dt_r_short_{}_{}_{}_{:07}.vtp".format(ID(E),int(N*NM_factor(E)),int(M*NM_factor(E)),time*time_factor(E))) for N,M in NMS] for E in EXPS])
for ix in range(1):#data_r2.shape[3]):
    print(data_e.shape)
    data_s = data_e[:,:,:,ix]
    sol_s = sol_r2[:,ix]
    plt.figure("Error  h_split -vs- h_reaction&h_diffusion combined. Component {}".format(ix))

    def err(data):
        return data - data[0,0,:]
    def err_true(data):
        return data - sol_s
    def norm_2(data):
        return np.linalg.norm(data, axis=2)
    def norm_oo(data):
        return np.linalg.norm(data, ord=np.inf, axis=2)
    def surp_2(data):
        return np.linalg.norm(data[1:,1:] - data[1:,:-1] - data[:-1,1:] + data[:-1,:-1], axis=2)
    def minmax(data):
        return np.min(data), np.max(data)

    err_s = err(data_s)
    err_s_2 = norm_2(err_s)
    err_s_oo = norm_oo(err_s)
    surp_s_2 = surp_2(data_s)

    cmap = matplotlib.cm.gray
    cmap.set_bad(color='red', alpha=0.3)
    cmap = matplotlib.cm.gray_r
    cmap.set_bad(color='red', alpha=0.3)


    plt.subplot(231)
    plt.imshow(nanan(np.log(np.abs(err_s_2))), cmap="gray_r", origin="lower")
    plt.title("||.||_2")
    plt.xlabel("NM factor")
    plt.ylabel("strang step")
    plt.xticks(np.arange(len(NMS)), NMS)
    plt.yticks(np.arange(len(EXPS)), EXPS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_e, (0,1)).T[::-1])
    plt.subplot(232)
    plt.title("||.||_oo")
    plt.imshow(nanan(np.log(np.abs(err_s_oo))), cmap="gray_r", origin="lower")
    plt.xticks(np.arange(len(NMS)), NMS)
    plt.yticks(np.arange(len(EXPS)), EXPS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_e, (0,1)).T[::-1])
    plt.subplot(233)
    plt.title("surplus")
    plt.imshow(nanan(np.log(np.abs(surp_s_2))), cmap="gray", origin="lower")
    plt.xticks(np.arange(len(NMS)-1), NMS)
    plt.yticks(np.arange(len(EXPS)-1), EXPS)
    plt.colorbar()










def SNM(S,N,M,base=None):
    if base is None: base=S
    if N%2 == 0 and M%4 == 0:
        return SNM(S - 1, N//2, M//2,  base)
    else:
        return S, N, M, time_factor(S, base)

def SNM_ID(S,N,M,base=None):
    S,N,M,tf = SNM(S, N, M, base)
    return ID(S), N, M, tf

NS=[256,128,64,32,16,8,4,2,1]
MS=[256,128,64,32,16,8,4,2]

data_dynamic_s = np.array([[try_load("out_neon_CN_magic_480/fibre_dt_r_short_{}_{}_{}_{:07}.vtp".format(*SNM_ID(1,N,M))) for N in NS] for M in MS])

for ix in range(1):#data_r2.shape[3]):
    data_s = data_dynamic_s[:,:,:,ix]
    sol_s  = data_dynamic_s[0,0,:,ix]
    plt.figure("Error h_reaction, h_diffusion, h_split = min possible. Component {}".format(ix))

    def err(data):
        return data - data[0,0,:]
    def err_true(data):
        return data - sol_s
    def norm_2(data):
        return np.linalg.norm(data, axis=2)
    def norm_oo(data):
        return np.linalg.norm(data, ord=np.inf, axis=2)
    def surp_2(data):
        return np.linalg.norm(data[1:,1:] - data[1:,:-1] - data[:-1,1:] + data[:-1,:-1], axis=2)
    def minmax(data):
        return np.min(data), np.max(data)

    err_s = err_true(data_s)
    err_s_2 = norm_2(err_s)
    err_s_oo = norm_oo(err_s)
    surp_s_2 = surp_2(data_s)

    mi,ma = minmax(np.array([minmax(nanan(np.log(np.abs(err_s_2)))), minmax(nanan(np.log(np.abs(err_d_2))))]))
    mima2 = {'vmin':mi, 'vmax':ma}
    print(mima2)

    plt.subplot(231)
    plt.imshow(nanan(np.log(np.abs(err_s_2))), cmap="gray_r", origin="lower", **mima2)
    plt.title("||.||_2")
    plt.xlabel("N (diffusion)")
    plt.ylabel("M (reaction)")
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_s, (0,0)).T[::-1])
    plt.plot(*opt_coarsen(data_s, (1,1)).T[::-1] + 0.1)
    plt.plot(*opt_refine(data_s).T[::-1] - 0.1)
    plt.subplot(232)
    plt.title("||.||_oo")
    plt.imshow(nanan(np.log(np.abs(err_s_oo))), cmap="gray_r", origin="lower")
    plt.xticks(np.arange(len(NS)), NS)
    plt.yticks(np.arange(len(MS)), MS)
    plt.colorbar()
    plt.plot(*opt_coarsen(data_s, (0,0), norm=np.inf).T[::-1])
    plt.plot(*opt_coarsen(data_s, (1,1), norm=np.inf).T[::-1] + 0.1)
    plt.plot(*opt_refine(data_s, norm=np.inf).T[::-1] - 0.1)
    plt.subplot(233)
    plt.title("surplus")
    plt.imshow(nanan(np.log(np.abs(surp_s_2))), cmap="gray", origin="lower")
    plt.xticks(np.arange(len(NS)-1), NS)
    plt.yticks(np.arange(len(MS)-1), MS)
    plt.colorbar()
# plt.show()




data_dynamic_s = np.array([[try_load("out_neon_CN_magic_480/fibre_dt_r_short_{}_{}_{}_{:07}.vtp".format(*SNM_ID(1,N,M))) for N in NS] for M in MS])
ns = np.array([[N for N in NS] for M in MS])
ms = np.array([[M for N in NS] for M in MS])
from mpl_toolkits.mplot3d import Axes3D
for ix in range(1):#data_r2.shape[3]):
    data_s = data_dynamic_s[:,:,:,ix]
    sol_s  = data_dynamic_s[0,0,:,ix]
    fig = plt.figure("Error h_reaction, h_diffusion, h_split = min possible. Component {}. 3D".format(ix))
    ax = fig.add_subplot(111, projection='3d')
    
    def err(data):
        return data - data[0,0,:]
    def err_true(data):
        return data - sol_s
    def norm_2(data):
        return np.linalg.norm(data, axis=2)
    def norm_oo(data):
        return np.linalg.norm(data, ord=np.inf, axis=2)
    def surp_2(data):
        return np.linalg.norm(data[1:,1:] - data[1:,:-1] - data[:-1,1:] + data[:-1,:-1], axis=2)
    def minmax(data):
        return np.min(data), np.max(data)

    err_s = err_true(data_s)
    err_s_2 = norm_2(err_s)
    err_s_oo = norm_oo(err_s)
    surp_s_2 = surp_2(data_s)
    
    ax.set_xlabel("Ns (diffusion)")
    ax.set_ylabel("Ms (raction)")
    ax.scatter(np.log2(ns),np.log2(ms),np.log2(err_s_2.reshape(-1)), label="||.||_2")
    ax.scatter(np.log2(ns),np.log2(ms),np.log2(err_s_oo.reshape(-1)), label="||.||_oo")
    ax.legend()
    
    # linear regression e = H c  =>  H^t H c = H^t e
    e = err_s_2.reshape(-1)
    # | hx_1^2 hy_1^2 hx_1 hy_1 |
    # | hx_2^2 hy_2^2 hx_2 hy_2 |
    # | hx_3^2 hy_3^2 hx_3 hy_3 |
    # |  ...    ...      ...    |
    H = 1.0/np.hstack([ns.reshape(-1,1)**2, ms.reshape(-1,1)**2, ns.reshape(-1,1)*ms.reshape(-1,1)])
    W = np.diag(np.sum(1/np.sqrt(H[:,:2]),axis=1))
    # c = np.linalg.solve(np.dot(H.T,H), np.dot(H.T,e))
    c = np.linalg.lstsq(H, e)[0]
    print('vvvvvvvvvv')
    print(c)
    c_log = np.linalg.solve(np.dot(H.T,np.dot(W,H)) + 0.0 * np.eye(H.shape[1]), np.dot(H.T,np.dot(W,e)))
    print(c_log)
    print('^^^^^^^^^^^^^')
    # ns = np.array([[N for N in 2.0**np.arange(-10,10)] for M in 2.0**np.arange(-10,10)])
    # ms = np.array([[M for N in 2.0**np.arange(-10,10)] for M in 2.0**np.arange(-10,10)])
    # ax.scatter(np.log2(ns),np.log2(ms),np.log2(c[0]*ns.reshape(-1)**-2.0 + c[1]*ms.reshape(-1)**-2.0 + c[2]*ns.reshape(-1)**-1.0*ms.reshape(-1)**-1.0), label="h1h1 h2h2 h1h2")
    # ax.scatter(np.log2(ns),np.log2(ms),np.log2(c[0]*ns.reshape(-1)**-2.0 + c[1]*ms.reshape(-1)**-2.0 +  0  *ns.reshape(-1)**-1.0*ms.reshape(-1)**-1.0), label="0*h1h2")
    # ax.scatter(np.log2(ns),np.log2(ms),c_log[0]*np.log2(ns.reshape(-1)**-2.0) + c_log[1]*np.log2(ms.reshape(-1)**-2.0) + c_log[2]*np.log2(ns.reshape(-1)**-1.0*ms.reshape(-1)**-1.0), label="logscale")
    # surf = ax.plot_surface(np.log2(ns),np.log2(ms),np.log2(c[0]*ns**-2.0 + c[1]*ms**-2.0 + c[2]*ns**-1.0*ms**-1.0), alpha=0.25, label="h1h1 h2h2 h1h2")
    # surf._facecolors2d=surf._facecolors3d
    # surf._edgecolors2d=surf._edgecolors3d
    surf = ax.plot_surface(np.log2(ns),np.log2(ms),np.log2(c_log[0]*ns**-2.0 + c_log[1]*ms**-2.0 + c_log[2]*ns**-1.0*ms**-1.0), alpha=0.25, label="weighted")
    # surf = ax.plot_surface(np.log2(ns),np.log2(ms),c_log[0]*np.log2((ns)**-2.0) + c_log[1]*np.log2((ms)**-2.0) + 0*c_log[2]*np.log2(ns**-1.0*2**ms**-1.0), alpha=0.25, label="logscale")
    surf._facecolors2d=surf._facecolors3d
    surf._edgecolors2d=surf._edgecolors3d
    
    ax.legend()
    
    
# plt.show()

xs = np.random.random(10)
np.sort(xs)
y = 5*xs*xs + 0.1*(np.random.random(10) - .5)
def kernel(xi,xj,eps=10, a = 0, q = 2):
    # return np.minimum(xi,xj) - xi * xj
    # return np.exp(-eps*np.abs(xi-xj))
    # return np.exp(-eps*np.abs(xi-xj)**2)
    return 1/np.maximum(xi,xj)**2 * np.minimum(xi,xj)**2
xi,xj = np.meshgrid(xs,xs)
K = kernel(xi,xj)
ai = np.linalg.solve(K, y)
s_xs = np.hstack([np.linspace(0,1,50)]);#, np.geomspace(1e-4,1,50)])
np.sort(s_xs)
s_ys = kernel(*np.meshgrid(xs, s_xs))
s = np.einsum("ij,j->i", s_ys, ai)
plt.figure("kernel")
plt.plot(s_xs, s)
plt.plot(s_xs, kernel(0.2,s_xs))
# plt.loglog(s_xs, kernel(0.2,s_xs))
# plt.loglog(s_xs, s)
plt.scatter(xs, y)

print(ai)
plt.show()
