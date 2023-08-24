 # https://github.com/t3rodrig/LammpsTutorials/blob/master/liquid-vapor/in.spce.txt
variable       Nsteps equal 500000    # 2 * 5e5 = 1e6 fs (1 ns). total number of steps in the simulation
variable       Nevery equal 100    # use input values every this many timesteps
variable       Nfreq equal ${Nsteps}/100    # calculate averages every this many timesteps
variable       Nrepeat equal 0.5*v_Nfreq/${Nevery}    # calculate averages every this many timesteps
variable       Text equal 300.0        # reference temperature, in K

processors     2 2 *
units          real
atom_style     full
boundary       p p p 

read_data      "data.long_box.interfaces"
include        "tip4p2005.ff"

reset_timestep 0

fix            fxshake all shake 0.001 20 0 b 1 a 1

timestep       0.5
fix            fxnvt all nvt temp ${Text} ${Text} $(150.0*dt)

thermo         1000
thermo_style   custom step dt cpu time temp press density pe ke
thermo_modify  flush yes
neigh_modify   every 1 delay 0 check yes

variable       xPress equal c_thermo_press[1] 
variable       yPress equal c_thermo_press[2] 
variable       zPress equal c_thermo_press[3]
variable       A_in_m equal 1e-10 # Angstrom in meter
variable       atm_in_Pa equal 101325 # note: 1 Pa = 1 N/m^2
variable       N_in_mN equal 1e3 # Newton in milliNewton
variable       convFac equal ${A_in_m}*${atm_in_Pa}*${N_in_mN}
variable       st equal 0.5*lz*(v_zPress-0.5*(v_xPress+v_yPress))*${convFac} # in mN/m
fix            st all ave/time 50 1 50 v_st file "out.surfTension.txt"
fix            strunning all ave/time 1 10 50 v_st ave running file "out.surfTension.running.txt"

variable       Dz equal 0.5 # Angstrom
compute        CC1 all chunk/atom bin/1d z center ${Dz} units box
fix            massDens all ave/chunk ${Nevery} ${Nrepeat} ${Nfreq} CC1 density/mass file "out.massDens.zProfile.txt"

# Computes the symmetric per-atom stress tensor
compute        CC2 all chunk/atom bin/1d z center ${Dz} units box
compute        stpa all stress/atom NULL
fix            stress all ave/chunk ${Nevery} ${Nrepeat} ${Nfreq} CC2 c_stpa[1] c_stpa[2] c_stpa[3] c_stpa[4] c_stpa[5] c_stpa[6] file "out.stress.zProfile.txt"

run            ${Nsteps}
write_data     "final.liquid-vapor.data" nocoeff
