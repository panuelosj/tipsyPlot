#%matplotlib inline
import numpy as np
import scipy.ndimage.filters as spfilt
import scipy.signal as spsignal
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import matplotlib.pyplot as plt
import sys, os, glob, pickle

 ######  ##          ###     ######   ######  ########  ######
##    ## ##         ## ##   ##    ## ##    ## ##       ##    ##
##       ##        ##   ##  ##       ##       ##       ##
##       ##       ##     ##  ######   ######  ######    ######
##       ##       #########       ##       ## ##             ##
##    ## ##       ##     ## ##    ## ##    ## ##       ##    ##
 ######  ######## ##     ##  ######   ######  ########  ######

class VariableProperties:
    ''' Holds properties for variables to be plotted '''

    def __init__(self, subplot=0, scatter=False, label=r'$NA$', title='NA'):
        self.subplot = subplot
        self.scatter = scatter
        self.label = label
        self.title = title

class Discontinuities:
    ''' Variable values over the discontinuous regions '''

    def __init__(self, r1=0, r2=0, r3=0, r4=0, r5=0):
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.r4 = r4
        self.r5 = r5

########     ###    ########     ###    ##     ##  ######
##     ##   ## ##   ##     ##   ## ##   ###   ### ##    ##
##     ##  ##   ##  ##     ##  ##   ##  #### #### ##
########  ##     ## ########  ##     ## ## ### ##  ######
##        ######### ##   ##   ######### ##     ##       ##
##        ##     ## ##    ##  ##     ## ##     ## ##    ##
##        ##     ## ##     ## ##     ## ##     ##  ######

# notes: remember to change the profile!!!
famname = 'rcrShock';  generictitle = 'RCR Shocktube\n'
nsteps = 100;           interval = 5
minxaxis = 0;           maxxaxis = 14
unitdistance = '1e3 kpc';   unitmass = 'Msol';  unittime = '1.5e10 Myr'
dMsolUnit = 1;          dKpcUnit = 1e3
dKmPerSecUnit = 6.55865e-05
dKmUnit = 3.085680859e19
dSecUnit = 4.70475e23;   dDelta = 0.01;
kmPerKpc = 3.241e-17;
dPressureUnit = (dKmPerSecUnit*dKmPerSecUnit)
xshift = 7
densityplotwidth = 28;   densityvmin = 0.001;  densityvmax = 0.125
profilebins = 200;
gamma = 1.4;
# diagnostics variables
threshold = 0.1

# format: 'name':{subplot(int), scatter(boolean)?, label name(str) }
variables = {'rho':VariableProperties(1, True, r'$\rho$', 'Density'),
             'p':VariableProperties(2, True, r'$P$', 'Pressure'),
             'u':VariableProperties(0, True, r'$u$', 'Internal Energy'),
             'velx':VariableProperties(4, True, r'$v_{x}$', 'Velocity (x component)'),
             'cSound':VariableProperties(0, True, r'$\sqrt{\frac{\gamma{} P}{\rho}}$', 'Sound Speed'),
             'entropy':VariableProperties(3, True, r'$\frac{P}{\rho{}^{\gamma{}}}$', 'Entropy')}





##     ## ########    ###    ########  ######## ########
##     ## ##         ## ##   ##     ## ##       ##     ##
##     ## ##        ##   ##  ##     ## ##       ##     ##
######### ######   ##     ## ##     ## ######   ########
##     ## ##       ######### ##     ## ##       ##   ##
##     ## ##       ##     ## ##     ## ##       ##    ##
##     ## ######## ##     ## ########  ######## ##     ##

# #Define stuffs
MIN_PEAK_HEIGHT = 0.05

# defined arrays
@pynbody.derived_array
def xaxis(sim):
    return sim['x'].in_units(unitdistance) + xshift
@pynbody.derived_array
def velx(sim):
    return sim['vel'][:,0]
@pynbody.derived_array
def cSound(sim):
    return np.sqrt((sim['p']*gamma)/(sim['rho']))
@pynbody.derived_array
def entropy(sim):
    return sim['p']/(sim['rho']**gamma)





######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######

def createPlot(varname, plotclear, scatter, title, savefig):
    '''
        Plots variables against the profile set for the simulation

        Arguments:
        varname     --  string
                        name of variable to be plotted
                        must correspond to the name assigned in-code (ie, density is 'rho')
        plotclear   --  boolean
                        if True, clears the window before plotting
                        generally True for plotting an individual variable
                            and False for plotting several variables in subplots
        scatter     --  boolean
                        if True, plots individual particle data on the same window
                            in addition to the profile's averaged data
        title       --  boolean
                        if True, plots a title above the plot based on the
                            title string assigned in the 'variables' dictionary
        savefig     --  boolean
                        if True, saves the finished figure to file with filename
                            [varname].[simulation].[simulationtime].png
    '''
    if plotclear == True:
        plt.clf()
    plt.axis([minxaxis, maxxaxis, globals()['min'+varname], globals()['max'+varname]])
    plt.plot(p['rbins'], p[varname])
    plt.xlabel('x [' + str(p['rbins'].units) + ']')
    plt.ylabel(variables[varname].label + ' [' + str(p[varname].units) + ']')

    # plot analytic solution
    plt.plot(arrx, globals()['arr' + varname])
    ## plot simulation discontinuities
    for i in range(0, len(simxgrad)):
        plt.axvline(x=simxgrad[i], color='r')
    for i in range(0, len(simxgrad2)):
        plt.axvline(x=simxgrad2[i], color='b')


    if scatter == True:
        plt.scatter(s['xaxis'], s[varname], marker='.')
    if title == True:
        plt.title(generictitle + variables[varname].title + ' (t=' + str(t).zfill(5) + ')')
    if savefig == True:
        plt.savefig('./test/' + varname + '/' + varname + '.' + simname + '.png',
                    dpi=100, bbox_inches='tight')

def updateMinMax(varname):
    '''
        A simple compare and update when looking for variable max and min values over time

        Arguments:
        varname     --  string
                        name of variable to be updated
                        must correspond to the name assigned in-code (ie, density is 'rho')
    '''
    if max(p[varname]) > globals()['max'+varname]:
        globals()['max'+varname] = max(p[varname])
    if min(p[varname]) < globals()['min'+varname]:
        globals()['min'+varname] = min(p[varname])







   ###     ######  ######## ##     ##    ###    ##           ######   #######  ########  ########
  ## ##   ##    ##    ##    ##     ##   ## ##   ##          ##    ## ##     ## ##     ## ##
 ##   ##  ##          ##    ##     ##  ##   ##  ##          ##       ##     ## ##     ## ##
##     ## ##          ##    ##     ## ##     ## ##          ##       ##     ## ##     ## ######
######### ##          ##    ##     ## ######### ##          ##       ##     ## ##     ## ##
##     ## ##    ##    ##    ##     ## ##     ## ##          ##    ## ##     ## ##     ## ##
##     ##  ######     ##     #######  ##     ## ########     ######   #######  ########  ########

nout = int(nsteps/interval)         # number of timesteps with simulation output
diag = np.zeros([nout,7])           # diagnostic values
#xTrue = np.zeros([nout,5])          # discontinuity positions based on the analytic solution
#xSim = np.zeros([nout,9])           # discontinuity positions from simulation
success = True

######## #### ##    ## ########     ########   #######  ##     ## ##    ## ########   ######
##        ##  ###   ## ##     ##    ##     ## ##     ## ##     ## ###   ## ##     ## ##    ##
##        ##  ####  ## ##     ##    ##     ## ##     ## ##     ## ####  ## ##     ## ##
######    ##  ## ## ## ##     ##    ########  ##     ## ##     ## ## ## ## ##     ##  ######
##        ##  ##  #### ##     ##    ##     ## ##     ## ##     ## ##  #### ##     ##       ##
##        ##  ##   ### ##     ##    ##     ## ##     ## ##     ## ##   ### ##     ## ##    ##
##       #### ##    ## ########     ########   #######   #######  ##    ## ########   ######
# find bounds   --  loop through the simulation (making profiles for each timestep) to find max and min values
for i in range(0,nout):
    t = (i+1)*interval
    simname = famname + '.' + str(t).zfill(5)
    s = pynbody.load(simname)
    s.physical_units(distance=unitdistance, mass=unitmass)
    p = profile.Profile(s, min=minxaxis, max=maxxaxis,
                        ndim=2, nbins=profilebins,
                        calc_x = lambda x: s.g['xaxis'])
    #p = profile.Profile(s)
    if i==0:
        for key in variables.keys():
            globals()['max'+key] = max(p[key]); globals()['min'+key] = min(p[key])
    else:
        for key in variables.keys():
            updateMinMax(key)

    if t%10 == 0:
        print('done t='+str(t).zfill(5))

for key in variables.keys():
    if globals()['max'+key] > 0:
        globals()['max'+key] *= 1.2
    else:
        globals()['max'+key] *= 0.8
    if globals()['min'+key] < 0:
        globals()['min'+key] *= 1.2
    else:
        globals()['min'+key] *= 0.8


   ###    ##    ##    ###    ##       ##    ## ######## ####  ######
  ## ##   ###   ##   ## ##   ##        ##  ##     ##     ##  ##    ##
 ##   ##  ####  ##  ##   ##  ##         ####      ##     ##  ##
##     ## ## ## ## ##     ## ##          ##       ##     ##  ##
######### ##  #### ######### ##          ##       ##     ##  ##
##     ## ##   ### ##     ## ##          ##       ##     ##  ##    ##
##     ## ##    ## ##     ## ########    ##       ##    ####  ######

 ######   #######  ##       ##     ## ######## ####  #######  ##    ##
##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ###   ##
##       ##     ## ##       ##     ##    ##     ##  ##     ## ####  ##
 ######  ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ##
      ## ##     ## ##       ##     ##    ##     ##  ##     ## ##  ####
##    ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ###
 ######   #######  ########  #######     ##    ####  #######  ##    ##

# before rarefaction, b/w rf tail and contact, b/w contact and rf tail, before rarefaction
truep = [0.083457531, 0.0, 0.0, 0.083457531]
truerho = [0.125, 0.0, 0.0, 0.125]
truevelx = [-5.0, 0.0, 0.0, 5.0]
truec = [0.0, 0.0, 0.0, 0.0]
trueentropy = [0.0, 0.0, 0.0, 0.0]
for i in range(0,len(truep)):
    try:
        truec[i] = np.sqrt(truep[i]*gamma/truerho[i])
    except ZeroDivisionError:
        truec[i] = 0.0

    try:
        trueentropy[i] = truep[i]/(truerho[i]**gamma)
    except ZeroDivisionError:
        trueentropy[i] = 0.0
trueu = [0.0, 0.0, 0.0, 0.0]

# before rarefaction, on rarefaction front, on rarefaction tail, on contact left, on contact right, on rarefaction tail, on rarefaction front, before rarefaction
arrp = [truep[0], truep[0], truep[1], truep[1], truep[2], truep[2], truep[3], truep[3]]
arrrho = [truerho[0], truerho[0], truerho[1], truerho[1], truerho[2], truerho[2], truerho[3], truerho[3]]
arrvelx = [truevelx[0], truevelx[0], truevelx[1], truevelx[1], truevelx[2], truevelx[2], truevelx[3], truevelx[3]]
arrcSound = [truec[0], truec[0], truec[1], truec[1], truec[2], truec[2], truec[3], truec[3]]
arrentropy = [trueentropy[0], trueentropy[0], trueentropy[1], trueentropy[1], trueentropy[2], trueentropy[2], trueentropy[3], trueentropy[3]]
arru = [trueu[0], trueu[0], trueu[1], trueu[1], trueu[2], trueu[2], trueu[3], trueu[3]]

# rarefaction front, rarefaction tail, contact, contact, rarefaction tail, rarefaction front
wavec = [-5.966811, -0.165942834, 0, 0, 0.165942834, 5.966811]

shockup = 4; shockdown = 3; rarefactionTail = 2; rarefactionFront = 1; contact = 0
wavestatus = [rarefactionFront, rarefactionTail, contact, contact, rarefactionTail, rarefactionFront]

arrp = [x*dPressureUnit for x in arrp]
arrvelx = [x*dKmPerSecUnit for x in arrvelx]
arrcSound = [x*dKmPerSecUnit for x in arrcSound]
arrentropy = [x*(dKmPerSecUnit**2) for x in arrentropy]
wavec = [x*dKmPerSecUnit*dSecUnit*dDelta*kmPerKpc*1e-3 for x in wavec]

wavec[1] = (arrvelx[1]+((2.0*arrcSound[1])/(0.4)))*dSecUnit*dDelta*kmPerKpc*1e-3
wavec[4] = (arrvelx[6]-((2.0*arrcSound[6])/(0.4)))*dSecUnit*dDelta*kmPerKpc*1e-3

print(wavec)


########  ##        #######  ########
##     ## ##       ##     ##    ##
##     ## ##       ##     ##    ##
########  ##       ##     ##    ##
##        ##       ##     ##    ##
##        ##       ##     ##    ##
##        ########  #######     ##

print('found bounds, creating plots now')
plt.figure(1)
plt.figure(2, figsize=(9,7))
plt.figure(3, figsize=(14,12))
plt.subplots(2, 2)
for i in range(0,nout):
    t = (i+1)*interval
    simname = famname + '.' + str(t).zfill(5)

    # import sim and make profile
    s = pynbody.load(simname)
    s.physical_units(distance=unitdistance, mass=unitmass)
    p = profile.Profile(s, ndim=2, nbins=profilebins, calc_x = lambda x: s.g['xaxis'])

    # analytic solution x's
    arrx = [minxaxis, wavec[0]*t + xshift, wavec[1]*t + xshift, wavec[2]*t + xshift, wavec[3]*t + xshift, wavec[4]*t + xshift, wavec[5]*t + xshift, maxxaxis]




    # find x's in simulation using first derivatives
    grad_rho = np.gradient(spfilt.gaussian_filter1d(p['rho'], 2.5), p['dr'][0])
    grad2_rho = np.gradient(spfilt.gaussian_filter1d(np.gradient(p['rho'], p['dr'][0]), 2.5), p['dr'][0])

    simxgrad = []
    simxgrad2 = []
    if t > 50:
        i = 0
        while i < len(wavestatus):
            if wavestatus[i] == contact:
                pass
            elif wavestatus[i] == rarefactionFront:
                # rarefaction front is where the profile is the most concave down,
                # therefore need to look for the second derivative's closest local minimum
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.less, order=2))
                simxgrad2.append(p['x'][abs(arrx[i+1]-xshift -p['x'])
                        == min(abs(arrx[i+1]-xshift -p['x'][arrextrema]))] +xshift)
            elif wavestatus[i] == rarefactionTail:
                # rarefaction tail is where the profile is the most concave up,
                # therefore need to look for the second derivative's closest local maximum
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.greater, order=2))
                simxgrad2.append(p['x'][abs(arrx[i+1]-xshift -p['x'])
                        == min(abs(arrx[i+1]-xshift -p['x'][arrextrema]))] +xshift)
            elif wavestatus[i] == shockup:
                # shock where the profile spikes up, therefore must look
                # for local max of the first derivative
                arrextrema = np.intersect1d(np.where(abs(grad_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad_rho, np.greater, order=2))
                simxgrad.append(p['x'][abs(arrx[i+1]-xshift -p['x'])
                        == min(abs(arrx[i+1]-xshift -p['x'][arrextrema]))] +xshift)

                # look for edges of the shock
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.greater, order=2))
                simxgrad2.append(p['x'][abs(simxgrad[len(simxgrad)-1]-xshift -p['x'])
                        == min(abs(simxgrad[len(simxgrad)-1]-xshift -p['x'][arrextrema]))] +xshift)
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.less, order=2))
                simxgrad2.append(p['x'][abs(simxgrad[len(simxgrad)-1]-xshift -p['x'])
                        == min(abs(simxgrad[len(simxgrad)-1]-xshift -p['x'][arrextrema]))] +xshift)
                i+=1
            elif wavestatus[i] == shockdown:
                # shock where the profile spikes down, therefore must look
                # for local min of the first derivative
                arrextrema = np.intersect1d(np.where(abs(grad_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad_rho, np.less, order=2))
                simxgrad.append(p['x'][abs(arrx[i+1]-xshift -p['x'])
                        == min(abs(arrx[i+1]-xshift -p['x'][arrextrema]))] +xshift)

                # look for edges of the shock
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.greater, order=2))
                simxgrad2.append(p['x'][abs(simxgrad[len(simxgrad)-1]-xshift -p['x'])
                        == min(abs(simxgrad[len(simxgrad)-1]-xshift -p['x'][arrextrema]))] +xshift)
                arrextrema = np.intersect1d(np.where(abs(grad2_rho)>MIN_PEAK_HEIGHT),
                        spsignal.argrelextrema(grad2_rho, np.less, order=2))
                simxgrad2.append(p['x'][abs(simxgrad[len(simxgrad)-1]-xshift -p['x'])
                        == min(abs(simxgrad[len(simxgrad)-1]-xshift -p['x'][arrextrema]))] +xshift)
                i+=1
            i+=1

    '''
    # place values into array
    xTrue[i] = [t, xAtrue, xBtrue, xCtrue, xDtrue]
    xSim[i] = [t, xAsim, xBsim, xCsim, xDsim, xCsimLower, xCsimUpper, xDsimLower, xDsimUpper]
    '''


    # plot gradient
    plt.figure(3)
    plt.clf()
    plt.suptitle(generictitle + '(t=' + str(t).zfill(5) + ')',
                 size='xx-large')
    plt.subplot(2, 2, 1)
    plt.xlim([minxaxis, maxxaxis])
    plt.autoscale(enable=True, axis='y')
    plt.plot(p['rbins'], grad_rho)
    plt.plot(p['rbins'], grad2_rho)

    ## plot simulation discontinuities
    for i in range(0, len(simxgrad)):
        plt.axvline(x=simxgrad[i], color='r')
    for i in range(0, len(simxgrad2)):
        plt.axvline(x=simxgrad2[i], color='b')
    plt.savefig('./test/gradients/gradients.' + simname + '.png', dpi=100, bbox_inches='tight')

    # export 2D slice density plot
    plt.figure(1)
    plt.clf()
    pp.image(s, width=densityplotwidth, vmin=densityvmin, vmax=densityvmax,
             units='1e-9 Msol kpc^-3')
    plt.suptitle(generictitle + 'Density Slice (t=' + str(t).zfill(5) + ')',
                 size='medium')
    plt.savefig('./test/fig/fig.' + simname + '.png', dpi=200, bbox_inches='tight')

    # Individual Plots
    plt.figure(2)
    for key in variables.keys():
        createPlot(key, True, variables[key].scatter, True, True)

    # create subplots
    plt.figure(3)
    plt.clf()
    plt.suptitle(generictitle + '(t=' + str(t).zfill(5) + ')',
                 size='xx-large')
    for key in variables.keys():
        if variables[key].subplot != 0:
            plt.subplot(2, 2, variables[key].subplot)
            #matplotlib.axes.ticklabel_format(style='sci', scilimits=(-2,2), axis='x')
            createPlot(key, False, variables[key].scatter, False, False)
    plt.savefig('./test/subplots/subplots.' + simname + '.png', dpi=100, bbox_inches='tight')



    if t%10 == 0:
        print('done t='+str(t).zfill(5) + '   p='+str(s['p']))

'''
# plot errors
diag[:,0] = xTrue[:,0]
# percent error of discontinuity positions
for i in range (1,5):
    diag[:,i] = ((xTrue[:,i]-xSim[:,i])/(xTrue[:,i]-xshift))
# discontinuity width as a percent of the simulation discontinuity position
diag[:,5] = ((xSim[:,6]-xSim[:,5])/(xSim[:,3]-xshift))
diag[:,6] = ((xSim[:,8]-xSim[:,7])/(xSim[:,4]-xshift))

# plot discontinuity error
plt.clf()
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,1], 'b-', label='error A')
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,2], 'g-', label='error B')
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,3], 'r-', label='error C')
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,4], 'c-', label='error D')
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,5], 'm-', label='width C')
plt.plot(diag[(np.ceil(200/interval)-1):,0], diag[(np.ceil(200/interval)-1):,6], 'y-', label='width D')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,1], 1), 'b--', label='smooth error A')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,2], 1), 'g--', label='smooth error B')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,3], 1), 'r--', label='smooth error C')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,4], 1), 'c--', label='smooth error D')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,5], 1), 'm--', label='smooth width C')
plt.plot(diag[(np.ceil(200/interval)-1):,0], spfilt.gaussian_filter1d(diag[(np.ceil(200/interval)-1):,6], 1), 'y--', label='smooth width D')
plt.legend(loc='upper right')
plt.xlabel('time [10 Ma]')
plt.ylabel('% error')
plt.title('Sod ShockTube Discontinuity Diagnostics')
plt.savefig('./diag/discontinuityerror.png', bbox_inches='tight')

plt.clf()
plt.plot(diag[(np.ceil(200/interval)-1):,0], np.gradient(spfilt.gaussian_filter1d(xSim[(np.ceil(200/interval)-1):,1], 1), interval), 'b-', label='ss A')
plt.plot(diag[(np.ceil(200/interval)-1):,0], np.gradient(spfilt.gaussian_filter1d(xSim[(np.ceil(200/interval)-1):,2], 1), interval), 'g-', label='ss B')
plt.plot(diag[(np.ceil(200/interval)-1):,0], np.gradient(spfilt.gaussian_filter1d(xSim[(np.ceil(200/interval)-1):,3], 1), interval), 'r-', label='ss C')
plt.plot(diag[(np.ceil(200/interval)-1):,0], np.gradient(spfilt.gaussian_filter1d(xSim[(np.ceil(200/interval)-1):,4], 1), interval), 'c-', label='ss D')

plt.plot([200,500], [-cs1,-cs1], 'b--', label='true A')
plt.plot([200,500], [cs3,cs3], 'g--', label='true B')
plt.plot([200,500], [cs4,cs4], 'r--', label='true C')
plt.plot([200,500], [cs5,cs5], 'c--', label='true D')
plt.legend(loc='upper right')
plt.xlabel('time [10 Ma]')
plt.ylabel('Speed')
plt.title('Sod ShockTube Wave Speed')
plt.savefig('./diag/wavespeed.png', bbox_inches='tight')
'''

# print animations
delay = 300/nout
os.system("convert -delay " + str(delay) + " ./test/fig/fig.*.png ./test/fig.animation.gif")
os.system("convert -delay " + str(delay) + " ./test/subplots/subplots.*.png ./test/subplots.animation.gif")
os.system("convert -delay " + str(delay) + " ./test/gradients/gradients.*.png ./test/gradients.animation.gif")
for key in variables.keys():
    os.system("convert -delay " + str(delay) + " ./test/"+key+"/"+key+".*.png ./test/"+key+".animation.gif")
