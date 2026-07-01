#
from abaqus import *
from abaqusConstants import *
from odbAccess import openOdb
import job, os, time, shutil, re, random, glob

# User configuration
SRC_INP = 'copper-1.inp'             # Source Abaqus INP file
JOB_NAME_BASE = 'copper-2'           # Target job name
ODB_NAME = 'copper-2.odb'            # Target ODB file
USER_SUB = 'E:\\Temp\\VUHARD-GE.f'   # Subroutine path
TOTAL_CYCLES = 10                    # Simulation iteration count

# Scan group file
def get_all_groups():
    node_files = glob.glob(os.path.join('NodeSetp-*.inp'))
    groups = {}
    for node_file in node_files:
        m = re.search(r'NodeSetp-(\d+)\.inp', os.path.basename(node_file))
        if m:
            idx = m.group(1)
            vf_file = os.path.join('VF-%s.inp' % idx)
            vel_file = os.path.join('Velocity-%s.inp' % idx)
            if os.path.exists(vf_file) and os.path.exists(vel_file):
                groups[idx] = {
                    'nodeset': node_file,
                    'vf': vf_file,
                    'velocity': vel_file
                }
    return groups

# Submit job
def submit_job(job_name):
    job_obj = mdb.jobs[job_name]
    job_obj.submit(consistencyChecking=OFF)
    job_obj.waitForCompletion()
    time.sleep(5)

# Extract ODB data
def extract_data(odb_file):
    odb = openOdb(path=odb_file)
    last_frame = odb.steps.values()[-1].frames[-1]

    # EVF output
    evf_field = last_frame.fieldOutputs['EVF_ASSEMBLY_EULERIAN_COPPER-1']
    with open('VF_output.inp', 'w') as f:
        for val in evf_field.values:
            if val.data>0:
                vf = 1.0 if val.data==1 else val.data
                f.write("Eulerian.%d, Eulerian.copper-1, %f\n" % (val.elementLabel, vf))

    # PEEQ output
    peeq_field = last_frame.fieldOutputs['PEEQ_ASSEMBLY_EULERIAN_COPPER-1']
    with open('PEEQ_output.inp','w') as f:
        for val in peeq_field.values:
            if abs(val.data)>1e-6:
                f.write("Eulerian.%d, %.6f\n" % (val.elementLabel, val.data))

    # Stress output
    stress_field = last_frame.fieldOutputs['S_ASSEMBLY_EULERIAN_COPPER-1']
    with open('Stress_output.inp','w') as f:
        for val in stress_field.values:
            if any(abs(c)>0.1 for c in val.data):
                data_str = ", ".join("%.6f" % x for x in val.data)
                f.write("Eulerian.%d, %s\n" % (val.elementLabel, data_str))

    # Temperature output
    temp_field = last_frame.fieldOutputs['NT11']
    with open('Temperature_output.inp','w') as f:
        for val in temp_field.values:
            if val.data>25:
                f.write("Eulerian.%d, %.6f\n" % (val.nodeLabel, val.data))

    # Velocity output
    vel_field = last_frame.fieldOutputs['V']
    with open('Velocity_output.inp','w') as f:
        for val in vel_field.values:
            v1,v2,v3 = val.data
            if abs(v1)>10000:
                f.write("Eulerian.%d,1,%.6f\n" % (val.nodeLabel, v1))
            if abs(v2)>10000:
                f.write("Eulerian.%d,2,%.6f\n" % (val.nodeLabel, v2))
            if abs(v3)>10000:
                f.write("Eulerian.%d,3,%.6f\n" % (val.nodeLabel, v3))
    odb.close()

# Modify INP file
def modify_inp_file(src_inp, group, first_cycle=False, last_cycle=False):
    nodeset_file = GROUPS[group]['nodeset']
    vf_file = GROUPS[group]['vf']
    velocity_file = GROUPS[group]['velocity']

    with open(src_inp,'r') as infile, open(JOB_NAME_BASE+'.inp','w') as outfile:
        skip_volume_fraction = False
        for line in infile:
            stripped = line.strip()
            
            # Modify the line following Dynamic Temperature-Displacement in the final iteration
            if last_cycle and '*Dynamic Temperature-displacement, Explicit' in line:
                outfile.write(line)
                next(infile)
                outfile.write(', 6e-07\n')
                continue

            # NodeSet
            if '*Nset, nset=NodeSet-p, instance=Eulerian' in line:
                outfile.write(line)
                outfile.write('*INCLUDE, INPUT=%s\n' % nodeset_file)
                continue
            # VF
            elif '*Initial Conditions, type=VOLUME FRACTION' in line:
                outfile.write(line)
                outfile.write('*INCLUDE, INPUT=%s\n' % vf_file)
                if not first_cycle:
                    outfile.write('*INCLUDE, INPUT=VF_output.inp\n')
                skip_volume_fraction = True
                continue
            elif skip_volume_fraction:
                if not first_cycle and re.match(r'^Eulerian\.(\d+), Eulerian\..*, ([\d.]+)', line.strip()):
                    continue
                else:
                    skip_volume_fraction=False
                    outfile.write(line)
                    continue
            # Temperature
            elif 'NodeSet-s, 25.' in line:
                outfile.write(line)
                if not first_cycle:
                    outfile.write('*INCLUDE, INPUT=Temperature_output.inp\n')
                continue
            # Velocity + Stress + PEEQ
            elif '*Initial Conditions, type=VELOCITY' in line:
                outfile.write(line)
                outfile.write('*INCLUDE, INPUT=%s\n' % velocity_file)
                if not first_cycle:
                    outfile.write('*INCLUDE, INPUT=Velocity_output.inp\n')
                    outfile.write('*INITIAL CONDITIONS, TYPE=STRESS\n')
                    outfile.write('*INCLUDE, INPUT=Stress_output.inp\n')
                    outfile.write('*INITIAL CONDITIONS, TYPE=HARDENING\n')
                    outfile.write('*INCLUDE, INPUT=PEEQ_output.inp\n')
                continue
            else:
                outfile.write(line)

# Main loop
def main_loop(total_cycles=TOTAL_CYCLES):
    global GROUPS
    GROUPS = get_all_groups()
    group_keys = list(GROUPS.keys())
    prev_seq = None
    sequence = group_keys[:]
    random.shuffle(sequence)

    for cycle in range(1, total_cycles+1):
        first_cycle = (cycle==1)
        last_cycle = (cycle == total_cycles)
        group = sequence[(cycle-1)%len(sequence)]

        # Generate copper-2.inp based on the source INP file
        modify_inp_file(SRC_INP, group, first_cycle=first_cycle, last_cycle=last_cycle)

        # Submit job
        job_current = JOB_NAME_BASE
        mdb.JobFromInputFile(name=job_current, inputFileName=JOB_NAME_BASE+'.inp', type=ANALYSIS,
                             memory=90, memoryUnits=PERCENTAGE, explicitPrecision=SINGLE,
                             nodalOutputPrecision=SINGLE, userSubroutine=USER_SUB,
                             resultsFormat=ODB, numDomains=16, activateLoadBalancing=False,
                             numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=16)
        submit_job(job_current)

        # Extract data (no extraction in the final iteration)
        if cycle < total_cycles:
            extract_data(ODB_NAME)
            # Modify ODB name
            new_odb_name = ODB_NAME.replace('.odb','(%d).odb' % cycle)
            shutil.move(ODB_NAME, new_odb_name)

if __name__=='__main__':
    main_loop()
