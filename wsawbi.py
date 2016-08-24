#!/usr/bin/python
# Copyright 2014-2016 Chris Sullivan
import readline, glob
import fortranformat as form
import numpy as np
import math
import pdb

def get(prompt, default):
    return raw_input("%s [%s] " % (prompt, default)) or default

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

def init_tab_complete():
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(complete)

def parse_dens(initialdens,finaldens):
    try:
        initial_dens_file = open(initialdens,"rb")
        final_dens_file = open(finaldens,"rb")
    except IOError:
        print "Error, file does not exist"; exit()

    for line in initial_dens_file:
        if '* ia,iz' in line:
            line = line.split()
            ai = int(line[3])
            zi = int(line[4])
            break
    for line in final_dens_file:
        if '* ia,iz' in line:
            line = line.split()
            af = int(line[3])
            zf = int(line[4])
            break
    print ai,zi," ",af,zf
    assert(ai==af)
    # determine which nucleon is valence in system
    if zi>zf:
        initial_nucleon = 'proton'
        final_nucleon = 'neutron'
    if zi<zf:
        initial_nucleon = 'neutron'
        final_nucleon = 'proton'
    if zi==zf:
        initial_nucleon = 'proton'
        final_nucleon = 'neutron'
        # doesn't matter initial == final
    valence_nucleon_initial = []
    sp_lines = False
    for i,line in enumerate(initial_dens_file):
        if sp_lines:
            if '-----' in line:
                sp_lines = False
                break
            valence_nucleon_initial.append(line.split())
        if initial_nucleon+' bound state results' in line:
            sp_lines = True
    valence_nucleon_final = []
    for i,line in enumerate(final_dens_file):
        if sp_lines:
            if '-----' in line:
                sp_lines = False
                break
            valence_nucleon_final.append(line.split())
        if final_nucleon+' bound state results' in line:
            sp_lines = True
    # print "initial valence nucleon: "
    # for line in valence_nucleon_initial:
    #     print line
    # print "final valence nucleon: "
    # for line in valence_nucleon_final:
    #     print line
    #return valence_nucleon_initial,valence_nucleon_final
    write_sp_levels(valence_nucleon_initial,ai,zi,initial_nucleon)
    write_sp_levels(valence_nucleon_final,af,zf,final_nucleon)


def write_sp_levels(levels_list,a,z,valence):
    corez = 0
    if valence == 'proton':
        corez = z - 1
    else:
        corez = z
    for i,lvl in enumerate(levels_list):
        if i == 1:
            continue # first line contains indices
        if len(lvl)==0:
            continue
        if len(lvl) != 9:
            print "Error in sp levels, bad level: ", lvl
            continue
        print i,lvl
        #line = form.FortranRecordWriter('(7F10.2,I2)')
        line = form.FortranRecordWriter('(F10.3,2F10.0,3F10.2,F10.1,I2)')
        line = line.write([a-1,corez,60.,0.65,1.25,1.25,7.0])
        file.write(line+'\n')
        line = form.FortranRecordWriter('(7F10.2,I2)')
        ebind = float(lvl[4])
        if ebind == -0.2:
            ebind = 1.0
        orbital_ang = lvl[2]
        if orbital_ang == 's':
            orbital_ang = 0
        elif orbital_ang == 'p':
            orbital_ang = 1
        elif orbital_ang == 'd':
            orbital_ang = 2
        elif orbital_ang == 'f':
            orbital_ang = 3
        elif orbital_ang == 'g':
            orbital_ang = 4
        elif orbital_ang == 'h':
            orbital_ang = 5
        else:
            print "error determinging orbital angular momentum: ", orbital_ang
            exit()
        ebind = abs(ebind)

        j = lvl[3]
        j = float(j[:-2])/float(j[len(j)-1])
        line = line.write([ebind, 1, orbital_ang, float(lvl[1])-1, 1 if valence == 'proton' else 0, j, 0.5])
        file.write(line+'\n')




if __name__=="__main__":
    #pdb.set_trace()
    init_tab_complete()
    #parse_obtd('t331dp150.obd')

    #FOLD = FOLDInstance()
    #FOLD.SetProjectile(13,5,1.5,-1,1.5)
    #FOLD.SetEjectile(13,4,0.5,+1,2.5)
    #FOLD.SetTarget(9,4,1.5,-1,0.5)
    #FOLD.SetResidual(9,5,1.5,-1,0.5)

    filename = get("Enter a name for the WSAW input file to be generated","wsaw.inp")
    file = open(filename,"wb")

    ################### Line 1 ###################

    line = form.FortranRecordWriter('(A10,A10,3I5)')
    line = line.write([0.1,20.,1,150,0])
    file.write(line+'\n')

    ################### Line 2 ###################
    output_file = \
        get("Enter a filename for the WSAW output file.\nNote that this filename is restricted to 8 characters or less","B13BE13")
    if len(output_file) > 8:
        print "Output WSAW filename is > 8 characters"
        exit()
    line = form.FortranRecordWriter('(A' + str(len(output_file))+ ')')
    line = line.write([output_file[0:8]])
    file.write(line+'\n')

    ################### Line 3-EOF ###################
    projdens = get("Enter path for dens projectile file from oxbash/nushell","13b-dens.dao")
    finaldens = get("Enter path for dens ejectile file from oxbash/nushell","13be-dens.dao")
    parse_dens(projdens,finaldens)

    line = form.FortranRecordWriter('(I2)')
    line = line.write([-1])
    file.write(line)
