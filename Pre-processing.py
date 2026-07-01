# Generate particles based on different distributions
from abaqus import *
from abaqusConstants import *
import visualization
import regionToolset
import mesh
import csv
import re

# =======================================
# USER CONFIGURATION - SET THESE PARAMETERS
MODEL_NAME = 'Model-1'  # Name of the model to work with
MERGED_PART_NAME = 'Particle'  # Name for the merged particle collection

# =======================================

def create_sphere_part(model, part_name, radius):
    """Create a spherical part with given name and radius"""
    s1 = model.ConstrainedSketch(name='__profile__', sheetSize=1000.0)
    s1.ConstructionLine(point1=(0.0, -0.1), point2=(0.0, 0.1))
    s1.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, radius),
                       point2=(0.0, -radius), direction=CLOCKWISE)
    s1.Line(point1=(0.0, radius), point2=(0.0, -radius))
    p = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    del model.sketches['__profile__']
    return p

# =======================================
def create_box_part(model, name, lx, ly, lz):
    s = model.ConstrainedSketch(name='__profile__', sheetSize=1000.0)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    p = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s, depth=lz)
    del model.sketches['__profile__']
    return p

# =======================================
def create_euler_part(model, name, lx, ly, lz):
    s = model.ConstrainedSketch(name='__profile__', sheetSize=1000.0)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    p = model.Part(name=name, dimensionality=THREE_D, type=EULERIAN)
    p.BaseSolidExtrude(sketch=s, depth=lz)
    del model.sketches['__profile__']
    return p

# =======================================
def read_sphere_data(csv_path):
    """Read sphere data from CSV file"""
    data_dict = {}
    with open(csv_path, 'rb') as file:
        reader = csv.reader(file)
        headers = reader.next()
        for header in headers:
            data_dict[header] = []
        for row in reader:
            for i, value in enumerate(row):
                data_dict[headers[i]].append(float(value))
    return data_dict
# =======================================
def process_particle_fraction(model, a, job, inp_file, i):
    pat = re.compile(r"Eulerian\.(\d+)", re.IGNORECASE)

    def get_vf_units(lines):
        vf, elems, flag = [], [], False
        for l in lines:
            if l.upper().startswith("*INITIAL CONDITIONS") and "VOLUME" in l.upper():
                flag = True
                continue
            if flag:
                if l.startswith("*"):
                    break
                vf.append(l)
                m = pat.search(l)
                if m: elems.append(int(m.group(1)))
        return vf, elems

    def get_nodes(elem_list, lines):
        if not elem_list:
            return []
        elem_set, nodes = set(elem_list), set()
        inside, target = False, False
        for l in lines:
            line = l.strip()
            if not line:
                continue
            if line.startswith("*"):
                inside = target = False
                if line.upper().startswith("*ELEMENT") and "EC3D8RT" in line.upper():
                    inside = target = True
                else:
                    inside = line.upper().startswith("*ELEMENT")
                continue
            if inside and target:
                parts = [p.strip() for p in line.split(",") if p.strip()]
                if len(parts) >= 2 and int(parts[0]) in elem_set:
                    nodes.update(map(int, parts[1:]))
        return sorted(nodes)

    def get_node_coords(lines, target_nodes):
        coords = {}
        target = set(target_nodes)
        inside = False
        for line in lines:
            s = line.strip()
            if s.upper().startswith("*NODE"):
                inside = True
                continue
            if inside and s.startswith("*"):
                break
            if inside:
                parts = [p.strip() for p in s.split(",") if p.strip()]
                if len(parts) >= 4:
                    nid = int(parts[0])
                    if nid in target:
                        coords[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                        if len(coords) == len(target):
                            break
        return coords

    def save(i, vf, nodes, node_coords):
        with open("VF-%02d.inp" % i, "w") as f:
            f.writelines(vf)
        with open("NodeSetp-%02d.inp" % i, "w") as f:
            for j in range(0, len(nodes), 16):
                f.write(", ".join(map(str, nodes[j:j + 16])) + "\n")
        with open("Velocity-%02d.inp" % i, "w") as f:
            for nid in sorted(node_coords):
                x, y, z = node_coords[nid]
                v = -615000
                f.write("Eulerian.%d, 3, %.6e\n" % (nid, v))

    # Read the Abaqus input file
    mdb.jobs[job].writeInput(consistencyChecking=OFF)
    lines = open(inp_file).readlines()
    vf, elems = get_vf_units(lines)
    nodes = get_nodes(elems, lines)
    node_coords = get_node_coords(lines, nodes)
    save(i, vf, nodes, node_coords)

# =======================================
# Configure the Eulerian domain, substrate geometry, and particle batches
def main(length_euler=0.35, width_euler=0.35, height_euler=0.6,
             height_substrate=0.05, height_particle=0.203, mesh_size=0.003,
             num_csv=5):
    
    model = mdb.models[MODEL_NAME]
    a = model.rootAssembly

    # Process CSV files in a loop
    for i in range(1, num_csv + 1):
        csv_file = "E:\\ABAQUS\\abaqus2019\\mutli particle create\\spheres_info-%02d.csv" % i

        # Delete particles from the previous iteration
        if i > 1:
            if 'Particle' in a.instances:
                del a.instances['Particle']

        data_dict = read_sphere_data(csv_file)
        x = data_dict['x']
        y = data_dict['y']
        z = data_dict['z']
        r = data_dict['r']
        num_particles = len(x)

        instances_info = []
        for j in range(num_particles):
            part_name = 'Part-{}'.format(j + 1)
            p = create_sphere_part(model, part_name, r[j])
            instances_info.append({
                'part': p,
                'position': (x[j], y[j], z[j])
            })

        temp_instances = []
        for j, info in enumerate(instances_info):
            instance_name = 'Temp-Instance-{}'.format(j + 1)
            a.Instance(name=instance_name, part=info['part'], dependent=ON)
            a.translate(instanceList=(instance_name,), vector=info['position'])
            temp_instances.append(instance_name)

        a.InstanceFromBooleanMerge(
            name=MERGED_PART_NAME,
            instances=[a.instances[name] for name in temp_instances],
            keepIntersections=ON,
            originalInstances=DELETE,
            domain=GEOMETRY
        )
        del a.instances['Particle-1']

        merged_part = model.parts[MERGED_PART_NAME]
        a.Instance(name="Particle", part=merged_part, dependent=ON)

        for j in range(num_particles):
            part_name = 'Part-{}'.format(j + 1)
            del model.parts[part_name]

        for key in a.features.keys():
            if 'Datum CSYS' in key:
                del a.features[key]

        # Additional initialization during the first iteration
        if i == 1:
            model = mdb.models[MODEL_NAME]

            euler_box = create_euler_part(model, 'Eulerian', length_euler, width_euler, height_euler)

            datum_plane = euler_box.DatumPlaneByPrincipalPlane(
                principalPlane=XYPLANE, offset=height_substrate)

            cells = euler_box.cells

            euler_box.PartitionCellByDatumPlane(
                datumPlane=euler_box.datums[datum_plane.id], cells=cells)

            matrix_box = create_box_part(model, 'Substrate', length_euler, width_euler, height_substrate)

            # Material assignment
            material = model.Material(name='Copper')
            model.materials['Copper'].Elastic(table=((128000.0, 0.34),))
            model.materials['Copper'].Conductivity(table=((385.0,),))
            model.materials['Copper'].Density(table=((8.96e-09,),))
            model.materials['Copper'].Plastic(hardening=USER, scaleStress=None, table=((0.0,),))
            model.materials['Copper'].InelasticHeatFraction(0.75)
            model.materials['Copper'].SpecificHeat(table=((384600000.0,),))
            model.EulerianSection(name='Section-1', data={'copper-1': 'Copper'})

            part = model.parts['Eulerian']
            region = part.Set(cells=part.cells[:], name='Set-1')
            part.SectionAssignment(region=region, sectionName='Section-1', offset=0.0,
                                   offsetType=MIDDLE_SURFACE, offsetField='',
                                   thicknessAssignment=FROM_SECTION)

            # Assembly
            a.Instance(name='Eulerian', part=euler_box, dependent=ON)
            a.Instance(name='Substrate', part=matrix_box, dependent=ON)
            a.translate(instanceList=['Eulerian'], vector=(0, 0, height_particle - height_euler))
            a.translate(instanceList=['Substrate'], vector=(0, 0, height_particle - height_euler))

            # Analysis step
            model.TempDisplacementDynamicsStep(name='Step-1', previous='Initial', timePeriod=3.25e-07,
                                               improvedDtMethod=ON)
            model.fieldOutputRequests['F-Output-1'].setValues(numIntervals=10, variables=(
                'S', 'PEEQ', 'V', 'NT', 'EVF'))
            del model.historyOutputRequests['H-Output-1']

            # Generate mesh
            p = model.parts['Eulerian']
            p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
            p.generateMesh()
            region = (p.cells,)
            p.setElementType(regions=region, elemTypes=(mesh.ElemType(elemCode=EC3D8RT, elemLibrary=EXPLICIT,
                                                                      secondOrderAccuracy=OFF,
                                                                      hourglassControl=DEFAULT),))

            # Boundary conditions
            f1 = a.instances['Eulerian'].faces
            faces1 = f1.getSequenceFromMask(mask=('[#4d2 ]',), )
            region = a.Set(faces=faces1, name='Set-1')
            model.VelocityBC(name='BC-1', createStepName='Initial',
                                             region=region, v1=0.0, v2=0.0, v3=0.0, vr1=UNSET, vr2=UNSET, vr3=UNSET,
                                             amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')
            f1 = a.instances['Eulerian'].faces
            faces1 = f1.getSequenceFromMask(mask=('[#24 ]',), )
            region = a.Set(faces=faces1, name='Set-2')
            model.VelocityBC(name='BC-2', createStepName='Initial',
                                             region=region, v1=0.0, v2=UNSET, v3=UNSET, vr1=UNSET, vr2=UNSET, vr3=UNSET,
                                             amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')
            f1 = a.instances['Eulerian'].faces
            faces1 = f1.getSequenceFromMask(mask=('[#108 ]',), )
            region = a.Set(faces=faces1, name='Set-3')
            model.VelocityBC(name='BC-3', createStepName='Initial',
                                             region=region, v1=UNSET, v2=0.0, v3=UNSET, vr1=UNSET, vr2=UNSET, vr3=UNSET,
                                             amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')

        # Execute in each iteration: create volume fraction, write inp, and process particle fraction
        a.regenerate()
        mdb.models[MODEL_NAME].rootAssembly.DiscreteFieldByVolumeFraction(
            name='DiscField-1', description='',
            eulerianInstance=a.instances['Eulerian'],
            referenceInstance=a.instances['Particle'])

        if i == 1:
            instList = (a.instances['Eulerian'],)
            rgn1 = a.instances['Eulerian'].sets['Set-1']
            fract1 = ("DiscField-1",)
            mdb.models[MODEL_NAME].MaterialAssignment(
                name='Predefined Field-1', instanceList=instList, useFields=True, fieldList=((rgn1, fract1),))

            # Modify the user subroutine path and the number of processors
            mdb.Job(name='copper-1', model=MODEL_NAME, description='', type=ANALYSIS,
                    memory=90, memoryUnits=PERCENTAGE,
                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
                    echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
                    userSubroutine='E:\\Temp\\VUHARD-GE.f',
                    resultsFormat=ODB, numDomains=16, numCpus=16)

        process_particle_fraction(model, a, 'copper-1', 'copper-1.inp', i)
        if 'DiscField-1' in a.features: del a.features['DiscField-1']
        
# =======================================
    a = model.rootAssembly
    a.regenerate()
    mdb.models['Model-1'].rootAssembly.DiscreteFieldByVolumeFraction(
        name='DiscField-1', description='',
        eulerianInstance=a.instances['Eulerian'],
        referenceInstance=a.instances['Substrate'])

    mdb.jobs['copper-1'].writeInput(consistencyChecking=OFF)

    pat = re.compile(r"Eulerian\.(\d+)", re.IGNORECASE)

    with open('copper-1.inp', 'r') as f:
        lines = f.readlines()

    vf_elems = []
    flag = False

    for l in lines:
        U = l.upper()
        if U.startswith("*INITIAL CONDITIONS") and "VOLUME" in U:
            flag = True
            continue

        if flag:
            if l.startswith("*"):  # End of the VF section
                break
            m = pat.search(l)
            if m:
                vf_elems.append(int(m.group(1)))

    vf_elem_set = set(vf_elems)
    nodes = []
    inside = False
    target = False

    for l in lines:
        line = l.strip()
        if not line:
            continue

        if line.startswith("*"):
            inside = target = False
            # Locate the EC3D8RT element block (same as above)
            if line.upper().startswith("*ELEMENT") and "EC3D8RT" in line.upper():
                inside = target = True
            else:
                inside = line.upper().startswith("*ELEMENT")
            continue

        if inside and target:
            parts = [p.strip() for p in line.split(",") if p.strip()]
            if len(parts) >= 2:
                eid = int(parts[0])
                if eid in vf_elem_set:
                    nodes.extend(map(int, parts[1:]))

    nodes = sorted(set(nodes))
    with open("NodeSets.inp", "w") as f:
        for i in range(0, len(nodes), 16):
            f.write(", ".join(map(str, nodes[i:i + 16])) + "\n")

# =======================================
    # Modify the input file after writing
    vf_line_pat = re.compile(r'^Eulerian\.\d+,\s*Eulerian\.copper-1,')

    with open('copper-1.inp', 'r') as f:
        lines = f.readlines()

    with open('copper-1.inp', 'w') as f:
        for i, line in enumerate(lines):
            stripped = line.strip()

            if '*Nset, nset=Set-1, instance=Eulerian' in line:
                f.write('*Nset, nset=NodeSet-p, instance=Eulerian\n')
                f.write('*Nset, nset=NodeSet-s, instance=Eulerian\n')
                f.write('*INCLUDE, INPUT=NodeSets.inp\n')

            f.write(line)

            # Set particle and substrate temperatures
            if vf_line_pat.match(stripped):
                if i + 1 == len(lines) or not vf_line_pat.match(lines[i + 1].strip()):
                    f.write('*Initial Conditions, type=TEMPERATURE\n')
                    f.write('NodeSet-p, 450.\n')
                    f.write('NodeSet-s, 25.\n')
                    f.write('*Initial Conditions, type=VELOCITY\n')

if __name__ == '__main__':
    main()
