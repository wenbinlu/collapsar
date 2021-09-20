import pylab as pl
import numpy as np
from math import sqrt, log10, pi, log, cos, floor
from scipy.interpolate import interp1d
from mplchange import *
import mesa_reader as mr
import sys

if len(sys.argv) > 2:
    fdir = sys.argv[1]     # directory for the files
    fname = sys.argv[2]    # file name supplied by command line
#    print('stellar model name: ' + fname)
    machine_run = True
    plt_history = False    # do not make plots
else:
    fdir = '/Users/wenbinlu/Documents/Research/collapsar/'   # manually supply the file directory
    fname = 'profile12'   # need to manually supply the file name
    machine_run = False
    plt_history = True


# adjustable parameters
s_PL = 0.5    # power-law index for the radial accretion rate profile in the advective regime
alp_ss = 0.03    # viscosity parameter
themin = 30*pi/180   # the minimum polar angle [rad] below which fallback is impeded by feedback
dt_over_tvis = 0.01       # time resolution factor

# some constants
c = 2.99792e10  # speed of light
G = 6.674e-8     # Newton's constant
rsun = 6.96e10         # solar radius
msun = 1.98847e33      # solar mass
R_unit = G*msun/c**2
J_unit = G*msun**2/c

prof = mr.MesaData(file_name=fdir+fname+'.data')
# print(prof.bulk_names)

rhodat = 10**np.flip(prof.logRho)       # density in each shell
rdat = np.flip(prof.radius)*rsun    # radius (right boundary)
Omgdat = np.flip(prof.omega)        # angular frequency of each radial shell

intp_lgrho = interp1d(np.log10(rdat), np.log10(rhodat), fill_value='extrapolate')
intp_lgomg = interp1d(np.log10(rdat), np.log10(Omgdat), fill_value='extrapolate')
# interpolate these profiles to a finer grid
Nr = 1000
rarr = np.logspace(log10(rdat[0]), log10(rdat[-1]), Nr)  # right boundary of shell
rhoarr = np.array([10**intp_lgrho(log10(r)) for r in rarr])
Omgarr = np.array([10**intp_lgomg(log10(r)) for r in rarr])


def risco_over_rg(a):   # ISCO radius
    z1 = 1 + (1-a*a)**(1./3) * ((1+a)**(1./3) + (1-a)**(1./3))
    z2 = (3*a*a + z1*z1)**0.5
    if a > 0:
        return 3+z2-((3-z1)*(3+z1+2*z2))**0.5
    return 3+z2+((3-z1)*(3+z1+2*z2))**0.5


def jisco_over_crg(a):      # specific AM at ISCO
    r = risco_over_rg(a)
    if a > 0:
        if a == 1:
            return (r**1.5 + r + r**0.5 - 1)/r**(3./4)/(r**0.5 + 2)**0.5
        return (r**2 - 2*a*r**0.5 + a**2)/r**(3./4)/(r**1.5 - 3*r**0.5 + 2*a)**0.5
    if a == -1:
        return (r**1.5 - r + r**0.5 + 1)/r**(3./4)/(r**0.5 - 2)**0.5
    return (r**2 + 2*a*r**0.5 + a**2)/r**(3./4)/(r**1.5 - 3*r**0.5 - 2*a)**0.5


Marr = np.empty(Nr, dtype=float)    # enclosed mass
ellarr = np.empty(Nr, dtype=float)  # specific angular momentum
Jarr = np.empty(Nr, dtype=float)    # enclosed angular momentum
tffarr = np.empty(Nr, dtype=float)  # free-fall time pi/2^1.5 * sqrt(r^3/GM)
rmid = np.empty(Nr, dtype=float)

M = 0.
J = 0.
for i in range(Nr):
    if i == 0:
        dM = rhoarr[i] * 4*pi/3 * rarr[i]**3
        dJ = 0.4 * Omgarr[i] * rarr[i]**2 * dM
    else:
        rhomid = (rhoarr[i] + rhoarr[i-1])/2
        dM = rhomid * 4*pi/3*(rarr[i]**3 - rarr[i-1]**3)
        omgmid = (Omgarr[i] + Omgarr[i-1])/2
        dJ = 4*pi * 2./3 * rhomid * omgmid * (rarr[i]**5 - rarr[i-1]**5)/5
    M += dM
    J += dJ
    Marr[i] = M
    Jarr[i] = J
    tffarr[i] = pi/2**1.5 * sqrt(rarr[i]**3/G/Marr[i])
    ellarr[i] = 2./3 * Omgarr[i] * rarr[i]**2

# need to find the time when accretion disk forms
Mbh0, abh0, Jbh0, i_disk = 0., 0., 0., 0
for i in range(Nr-1):
    if Marr[i] < msun:   # assume that the innermost 1 Msun always forms a BH
        continue
    Mbh = Marr[i]
    Rg = G*Mbh/c**2
    abh = c*Jarr[i]/(G*Mbh**2)
    ell_in_crg = ellarr[i+1]/c/Rg
    jisco_in_crg = jisco_over_crg(abh)
    if ell_in_crg > jisco_in_crg:
        i_disk = i+1
        Mbh0 = Marr[i]
        Jbh0 = Jarr[i]
        abh0 = abh
        break
    if i == Nr-2:
        if machine_run:
            print('%.3e\t%.3e\t%.3e\t%.3e\t%.3e' % (Mbh/msun, abh, 0, 0, 0))
        else:
            print('no disk forms for this star!')
            print('final Mbh=%.3f Msun, abh=%.3f' % (Mbh/msun, abh))
        exit()
# the i_disk-th shell starts forming a disk
tdisk = tffarr[i_disk]

# fallback rate profile
fb_frac = cos(themin)   # fraction of mass outside 30-deg polar cones
Mfbdot = np.diff(Marr)/np.diff(tffarr) * fb_frac
Jfbdot = np.diff(Jarr)/np.diff(tffarr) * (cos(themin) - 1./3 * (cos(themin))**3) * 3./2
tmid = np.array([(tffarr[i]+tffarr[i+1])/2 for i in range(Nr-1)])
intp_lgMfbdot = interp1d(tmid, np.log10(Mfbdot), fill_value='extrapolate')
intp_lgJfbdot = interp1d(tmid, np.log10(Jfbdot), fill_value='extrapolate')

tmax = tffarr[-1]*1.5
Ntgrid = 300        # interpolate lgMfbdot and lgJfbdot on a regular grid
tgrid = np.linspace(tdisk, tmax, Ntgrid)
dtgrid = tgrid[1] - tgrid[0]
lgMfbdotgrid = np.empty(Ntgrid, dtype=float)
lgJfbdotgrid = np.empty(Ntgrid, dtype=float)
for i in range(Ntgrid):
    t = tgrid[i]
    if t > tmid[-1]:
        lgMfbdotgrid[i] = -10.   # zero fallback rate
        lgJfbdotgrid[i] = -10.
    else:
        lgMfbdotgrid[i] = intp_lgMfbdot(t)
        lgJfbdotgrid[i] = intp_lgJfbdot(t)


def MJfbdot(t, tgrid, lgMJdot_grid):    # M/J fallback rate any time
    i_grid = min(Ntgrid-1, max(0, int(floor((t-tgrid[0])/dtgrid))))
    # t is usually between tgrid[i_grid] and tgrid[i_grid+1]
    slope = (lgMJdot_grid[i_grid+1] - lgMJdot_grid[i_grid])/(tgrid[i_grid+1] - tgrid[i_grid])
    lgMJdot = lgMJdot_grid[i_grid] + (t - tgrid[i_grid])*slope
    return 10**lgMJdot


tarr = []       # time
Mdarr = []      # disk mass
Rdarr = []      # outer disk radius
Riscoarr = []    # ISCO radius
Mbharr = []    # BH mass
Mfbdotarr = []  # mass fallback rate
Jfbdotarr = []  # AM fallback rate
Mbhdotarr = []   # BH mass gaining rate
Mdotaccarr = []   # outer disk accretion rate
Liscoarr = []    # accretion power near isco
Lnuiscoarr = []     # neutrino power near isco
Lwiscoarr = []      # wind power near isco
Lwarr = []    # total wind power for the entire disk

Mbh = Mbh0
Jbh = Jbh0
abh = c*Jbh/(G*Mbh**2)
Rg = G*Mbh/c**2
Risco = risco_over_rg(abh) * Rg
OmgKisco = sqrt(G*Mbh/Risco**3)
tvis_isco = 1/alp_ss/OmgKisco

# initialize the disk properties (unimportant for the total energetics)
Md0 = MJfbdot(tdisk, tgrid, lgMfbdotgrid)*tvis_isco*0.5
Rd0 = 1.03*Risco
Jd0 = sqrt(G*Mbh*Rd0) * Md0   # initial disk amgular momentum

Md = Md0
Jd = Jd0
t = tdisk
while t < tmax:
    abh = c*Jbh/(G*Mbh**2)
    Rg = G*Mbh/c**2
    Risco = risco_over_rg(abh) * Rg
    Rd = (Jd/Md)**2/(G*Mbh)   # Newtonian Keplerian rotation
    OmgK = sqrt(G*Mbh/Rd**3)
    tvis = 1/alp_ss/OmgK
    Mdotacc = Md/tvis
    Rt = max(Risco, min(Rd, (2*Rg/Rd**s_PL * Mdotacc/msun * 10**2.5)**(1./(1-s_PL))))
    Mbhdot = Mdotacc * (Rt/Rd)**s_PL
    Jwdot = 2*s_PL/(2*s_PL + 1) * sqrt(G*Mbh*Rd) * Mdotacc *\
        (1 - (Rt/Rd)**((2*s_PL+1)/2))
    Jbhdot = jisco_over_crg(abh) * Rg * c * Mbhdot

    Mfbdot = MJfbdot(t, tgrid, lgMfbdotgrid)
    Jfbdot = MJfbdot(t, tgrid, lgJfbdotgrid)

    Mddot = Mfbdot - Mdotacc
    Jddot = Jfbdot - Jwdot - Jbhdot
    Lw = 0.5*s_PL/(1-s_PL)*G*Mbh/Rd * Mdotacc * ((Rd/Rt)**(1-s_PL) - 1)
    eta_acc = 1 - sqrt(1 - 2*Rg/3./Risco)
    Lacc = Mbhdot*c**2*eta_acc
    if Rt > 1.001*Risco:
        Lnuisco = Lacc
        Lwisco = 0.
    else:
        Lnuisco = 0.
        Lwisco = Lacc

    tarr += [t]
    Mdarr += [Md/msun]
    Rdarr += [Rd]
    Riscoarr += [Risco]
    Mbharr += [Mbh/msun]
    Mfbdotarr += [Mfbdot/msun]
    Jfbdotarr += [Jfbdot/J_unit]
    Mbhdotarr += [Mbhdot/msun]
    Mdotaccarr += [Mdotacc/msun]
    Lwarr += [Lw]
    Liscoarr += [Lacc]
    Lnuiscoarr += [Lnuisco]
    Lwiscoarr += [Lwisco]
    dt = dt_over_tvis * tvis
    t += dt
    Md += Mddot * dt
    Jd += Jddot * dt
    Mbh += Mbhdot * dt
    Jbh += Jbhdot * dt

Nt = len(tarr)

Ewind = 0.      # total wind energy
Eacc = 0.        # total accretion energy
Eacc_ADAF = 0.   # accretion energy when isco is ADAF
for i in range(Nt-1):
    dt = tarr[i+1] - tarr[i]
    Ewind += Lwarr[i] * dt
    Eacc += Liscoarr[i] * dt
    Eacc_ADAF += Lwiscoarr[i] * dt

if machine_run:
    print('%.3e\t%.3e\t%.3e\t%.3e\t%.3e' % (Mbh/msun, abh, Ewind, Eacc, Eacc_ADAF))
else:
    print('disk formation time since core-collapse: %.3f s' % tdisk)
    print('total wind energy: Ew=%.3e erg' % Ewind)
    print('total accretion energy: %.3e erg' % Eacc)
    print('accretion energy when ISCO is advection-dominated (for jet formation): %.3e erg' % Eacc_ADAF)
    print('final Mbh=%.3f Msun, abh=%.3f' % (Mbh/msun, abh))

# the rest of the code makes some plots
if not plt_history:
    exit()

tarr = np.array(tarr)
# print('number of timesteps: %d' % Nt)
tplt = tarr - tdisk
fig, (ax3, ax1, ax2) = pl.subplots(3, 1, figsize=(13, 18), sharex='all')
leg_fs = 30

ax3.plot(tplt, Rdarr, 'r-', lw=3, label=r'$R_{\rm d}$', alpha=0.5)
ax3.plot(tplt, Riscoarr, 'b--', lw=3, label=r'$R_{\rm isco}$', alpha=0.5)
ax3.set_ylabel(r'$R\rm\,[cm]$', labelpad=-2)
ax3.set_yscale('log')
ax3.legend(loc='best', prop={'size': leg_fs}, fancybox=True, framealpha=0.3)

ax1.plot(tplt, Lwarr, 'r-', lw=3, label=r'$L_{\rm w,tot}$', alpha=0.5)
ax1.plot(tplt, Lnuiscoarr, 'b--', lw=3, label=r'$L_{\rm\nu, isco}$', alpha=0.5)
ax1.plot(tplt, Lwiscoarr, 'g:', lw=3, label=r'$L_{\rm w, isco}$', alpha=0.8)
ax1.set_ylim(1.73e-4*max(Liscoarr), 1.72*max(Liscoarr))
ax1.set_ylabel(r'$L\rm\,[erg\,s^{-1}]$', labelpad=-2)
ax1.set_yscale('log')
ax1.legend(loc='best', prop={'size': leg_fs}, fancybox=True, framealpha=0.3)

ax2.plot(tplt, Mfbdotarr, 'r-', lw=3, ms=3,
        label=r'$\dot{M}_{\rm fb}$', alpha=0.5)
ax2.plot(tplt, Mdotaccarr, 'b--', lw=3, ms=3,
        label=r'$\dot{M}_{\rm acc}$', alpha=0.5)
ax2.plot(tplt, Mbhdotarr, 'g:', lw=3, ms=3,
         label=r'$\dot{M}_{\rm BH}$', alpha=0.8)

ax2.set_yscale('log')
ax2.set_ylim(3.5e-6, 0.35)
ax2.set_ylabel(r'$\dot{M}\,[M_\odot\rm\,s^{-1}]$', labelpad=-2)
ax2.legend(loc='best', prop={'size': leg_fs}, fancybox=True, framealpha=0.3)

ax2.set_xscale('log')
ax2.set_xlim(2.5e-3, 550)
ax2.set_xlabel(r'$t-t_{\rm disk}\rm\,[s]$', labelpad=-2)

pl.tick_params(axis='both', which='both', direction='in', bottom=True,
               top=True, left=True, right=True)
pl.subplots_adjust(wspace=0, hspace=0)
pl.subplots_adjust(bottom=0.07, left=0.12, top=0.98, right=0.98)
fig.savefig(fdir + 'evolution_' + fname + '.png')
