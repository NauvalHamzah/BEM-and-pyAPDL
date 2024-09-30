import os
import numpy as np
from ansys.mapdl.core import launch_mapdl

#mapdl = launch_mapdl()

'''
new_dir = 'F:/Kuliah/Tesis/Program/ANSYS_APDL/NREL5MW/'
fname = new_dir+'Trial'
inp_command = 'aTransientCommand'

mapdl.cwd(new_dir)
mapdl.clear()
mapdl.resume('Trial','db')
#mapdl.input(inp_command,'inp')
mapdl.slashsolu()
mapdl.omega(0,0.1,0)
mapdl.acel(0,0,1)
mapdl.lssolve(1,1,1)
#mapdl.save("Trial",'db')
'''


def cylinder_batch(elemsize, plot=False):
    """ Report the maximum von Mises stress of a Cantilever supported cylinder"""
    # clear
    mapdl.finish()
    mapdl.clear()
    # cylinder parameters
    radius = 2
    h_tip = 2
    height = 20
    force = 200/radius
    pressure = force/(h_tip*2*np.pi*radius)

    mapdl.prep7()
    mapdl.et(1, 186)
    mapdl.et(2, 154)
    mapdl.r(1)
    mapdl.r(2)

    # Aluminum properties (or something)
    mapdl.mp('ex', 1, 10e6)
    mapdl.mp('nuxy', 1, 0.3)
    mapdl.mp('dens', 1, 0.1/386.1)
    mapdl.mp('dens', 2, 0)

    # Simple cylinder
    for i in range(4):
        mapdl.cylind(radius, '', '', height, 90*(i-1), 90*i)

    mapdl.nummrg('kp')

    # mesh cylinder
    mapdl.lsel('s', 'loc', 'x', 0)
    mapdl.lsel('r', 'loc', 'y', 0)
    mapdl.lsel('r', 'loc', 'z', 0, height - h_tip)

    # mapdl.lesize('all', elemsize*2)
    mapdl.mshape(0)
    mapdl.mshkey(1)
    mapdl.esize(elemsize)
    mapdl.allsel('all')
    mapdl.vsweep('ALL')
    mapdl.csys(1)
    mapdl.asel('s', 'loc', 'z', '', height - h_tip + 0.0001)
    mapdl.asel('r', 'loc', 'x', radius)
    mapdl.local(11, 1)
    mapdl.csys(0)
    mapdl.aatt(2, 2, 2, 11)
    mapdl.amesh('all')
    mapdl.finish()

    if plot:
        mapdl.view(1, 1, 1, 1)
        mapdl.eplot()

    # new solution
    mapdl.slashsolu()
    mapdl.antype('static', 'new')
    mapdl.eqslv('pcg', 1e-8)

    # Apply tangential pressure
    mapdl.esel('s', 'type', '', 2)
    mapdl.sfe('all', 2, 'pres', '', pressure)

    # Constrain bottom of cylinder/rod
    mapdl.asel('s', 'loc', 'z', 0)
    mapdl.nsla('s', 1)
    mapdl.run('d, all , all')
    mapdl.allsel()
    mapdl.psf('pres', '', 2)
    mapdl.pbc('u', 1)
    mapdl.run('solve')
    #mapdl.save('cylinder','db')
    mapdl.finish()

    # access results using MAPDL object
    result = mapdl.result

    # to access the results you could have run:
    # from ansys.mapdl import reader as pymapdl_reader
    # resultfile = os.path.join(mapdl.path, '%s.rst' % mapdl.jobname)
    # result = pymapdl_reader.read_binary(result file)
    # Get maximum von Mises stress at result 1
    # Index 0 as it's zero based indexing
    nodenum, stress = result.principal_nodal_stress(0)

    # von Mises stress is the last column
    # must be nanmax as the shell element stress is not recorded
    maxstress = np.nanmax(stress[:, -1])

    # return number of nodes and max stress
    return nodenum.size, maxstress
    
def resume (elemsize,dbfile,new_dir, plot=False):
    # clear
    mapdl.finish()
    mapdl.clear()
    
    mapdl.filname('Blade')
    # resume dbfile
    resume_command = 'resume,'+dbfile+',db'
    mapdl.run(resume_command)
    mapdl.run('/input,F:/Kuliah/Tesis/Program/ANSYS_APDL/Cylinder/input_apdl.txt')
    mapdl.run('/input,F:/Kuliah/Tesis/Program/ANSYS_APDL/Cylinder/cyl_out.txt')

    '''
    mapdl.run('/post1')
    mapdl.run('/output,F:/Kuliah/Tesis/Program/ANSYS_APDL/Cylinder/Blade_Result.txt')
    mapdl.run('PRNSOL,U,SUM')
    mapdl.run('/out')
    '''

    nodenum=1
    maxstress=2
    # return number of nodes and max stress
    return nodenum, maxstress


# initialize MAPDL
mapdl = launch_mapdl(override=True, loglevel='ERROR')
new_dir = 'F:/Kuliah/Tesis/Program/ANSYS_APDL/Cylinder/'
new_dir2 = 'C:/Temp'
new_db = new_dir+"/Blade"
new_db2 = new_dir2+"/Trial"

result_summ = []
elemsize=2
nnode, maxstress = resume(elemsize,new_db2,new_dir, plot=False)
result_summ.append([nnode, maxstress])
print('Element size %f: %6d nodes and maximum vom Mises stress %f'% (elemsize, nnode, maxstress))


'''
for elemsize in np.linspace(2, 0.15, 1):
    # run the batch and report the results
    nnode, maxstress = cylinder_batch(elemsize, plot=False)
    result_summ.append([nnode, maxstress])
    print('Element size %f: %6d nodes and maximum vom Mises stress %f'
    % (elemsize, nnode, maxstress))
'''
mapdl.save(new_db,'db')
# Exit MAPDL
mapdl.exit()

