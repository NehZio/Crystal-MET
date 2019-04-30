# Writing an INPUT file

There is no obligation to put any order in the keywords

## MANDATORY KEYWORDS

1.BATH bath_radius

1.PSEUDO first_shell_radius

1.OUTPUT output_file_name

1.PATTERN nAtom1 labelAtom1 nAtom2 labelAtom2 ...  

Not that, if you have more than once the Atom1, you will have to place it again at the end (for example, 5 Ir 8 Sr 6 O 5 Ir)

1.NPATTERN nPattern

1.LATTICE
  a cell_parameter_a
  b cell_parameter_b
  c cell_parameter_c
  alpha cell_angle_alpha
  beta cell_angle_beta
  gamma cell_angle_gamma

1. ATOMS label charge neighbours radius_of_the_coordination_sphere ...

1.SYMOP
   x,y,z
   ...
   POMYS
   
1.GENERATOR
  atom1 x_frac y_frac z_frac
  ...
  ROTARENEG
  
1.END_OF_INPUT


## OPTIONNAL KEYWORDS

1.TRANSLATE x y z

Crystal MET has some trouble finding specific pattern, sometimes forcing a translation can help finding it

1.CENTER  atom1 atom2 ...

Will chose the atom (or the center of the atoms) as the center of the sphere

1.X_AXIS atom1 atom2 ...

Will chose the atom (or the center of the atoms) as the x axis

1.Y_AXIS  atom1  atom2 ...

Will chose the atom (or the center of the atoms) as the y axis

1.Z_AXIS atom1  atom2 ...

Will chose the atom (or the center of the atoms) as the z axis

1.NOTINPP atoms not in the first shell

1.COLOR  

Shows the bath sorted with colors

1.NOCOLOR 

Shows the bath sorted by atoms

1.SEEFRAG  

Shows the fragment once found (recommended)

1.OPTIMIZATION  

Will use an optimization method to set the total charge to 0 by slightly changing the charges of the outer shell atoms

1.EVJEN 
Simili Evjen method to calculate the charge of the atom on the outer shell, note that OPTIMIZATION and EVJEN keywords can not be used together

