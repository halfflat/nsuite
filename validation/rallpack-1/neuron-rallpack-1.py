#!/usr/bin/env python

from math import pi
import sys
import os

from neuron import h
import nsuite.stdarg as stdarg
import xarray

ra =     1.0   # axial resistivity [Ω m]
rm =     4.0   # membrane resistivity [Ω m²]
cm =    0.01   # memrane specific capacitance [F/m²]
Erev =   -65   # reversal potential [mV]

diam =    1.0  # cable diameter [µm]
length = 1000  # cable length [µm]
iinj =    0.1  # current injection [nA]
n =      1000  # compartments
x = 333.33333  # measurement point [µm]

dt = 0.0025

sample_dt = 0.001
tend = 0.25

# One recognized tag: 'firstorder'

output, tags, params = stdarg.parse_run_stdarg(tagset=['firstorder'])
for v in ['dt', 'n']:
    if v in params: globals()[v] = params[v]

# TODO: resolve routines below with exisiting neuron support code.

def hoc_execute_quiet(arg):
    with open(os.devnull, 'wb') as null:
        fd = sys.stdout.fileno()
        keep = os.dup(fd)
        sys.stdout.flush()
        os.dup2(null.fileno(), fd)
        h(arg)
        sys.stdout.flush()
        os.dup2(keep, fd)

def hoc_setup():
    hoc_execute_quiet('load_file("stdrun.hoc")')

# Make model

hoc_setup()

cable = h.Section(name='cable')
cable.diam = diam
cable.L= length
cable.cm = 100*cm       # [µF/cm² = 0.01 F/m²]
cable.Ra = 100*ra       # [Ω cm = 0.01 Ω m]
cable.nseg = 1000

cable.insert('pas')
cable.g_pas = 10000/rm  # [S/cm² = 0.0001 S/m²]
cable.e_pas = Erev

h.v_init = Erev

# Run model

v = h.Vector()
v.record(cable(1.0/3.0)._ref_v, sample_dt)
t = h.Vector()
t.record(h._ref_t, sample_dt)

h.dt = dt
h.steps_per_ms = 1/dt # or else NEURON might noisily fudge dt
if 'firstorder' in tags:
    h.secondorder = 0
else:
    h.secondorder = 2
h.tstop = tend
h.run()

# Collect and save data

out = xarray.Dataset({'voltage': (['time'], list(v))}, coords={'time': list(t)})
out['n'] = n
out['dt'] = dt
out.to_netcdf(output)

