load ref_target.pdb, ref
hide everything, ref
show cartoon, ref
color gray80, ref
bg_color white

# pocket center
pseudoatom p0, pos=[38.379,34.977,36.580]
show spheres, p0
set sphere_scale, 0.5, p0
color yellow, p0

# arrows for 10 trials
set dash_width, 4
set dash_gap, 0
set label_size, 20
pseudoatom pEnd01, pos=[19.754,20.107,33.749]
distance dir01, p0, pEnd01
color blue, dir01
show spheres, pEnd01
set sphere_scale, 0.25, pEnd01
color blue, pEnd01
label pEnd01, "1"
pseudoatom pEnd02, pos=[16.290,26.396,32.781]
distance dir02, p0, pEnd02
color blue, dir02
show spheres, pEnd02
set sphere_scale, 0.25, pEnd02
color blue, pEnd02
label pEnd02, "2"
pseudoatom pEnd03, pos=[26.324,29.274,16.626]
distance dir03, p0, pEnd03
color blue, dir03
show spheres, pEnd03
set sphere_scale, 0.25, pEnd03
color blue, pEnd03
label pEnd03, "3"
pseudoatom pEnd04, pos=[17.114,24.222,39.437]
distance dir04, p0, pEnd04
color blue, dir04
show spheres, pEnd04
set sphere_scale, 0.25, pEnd04
color blue, pEnd04
label pEnd04, "4"
pseudoatom pEnd05, pos=[58.756,24.946,44.335]
distance dir05, p0, pEnd05
color red, dir05
show spheres, pEnd05
set sphere_scale, 0.25, pEnd05
color red, pEnd05
label pEnd05, "5"
pseudoatom pEnd06, pos=[19.339,26.287,24.834]
distance dir06, p0, pEnd06
color blue, dir06
show spheres, pEnd06
set sphere_scale, 0.25, pEnd06
color blue, pEnd06
label pEnd06, "6"
pseudoatom pEnd07, pos=[19.685,22.349,28.389]
distance dir07, p0, pEnd07
color blue, dir07
show spheres, pEnd07
set sphere_scale, 0.25, pEnd07
color blue, pEnd07
label pEnd07, "7"
pseudoatom pEnd08, pos=[16.909,29.635,27.278]
distance dir08, p0, pEnd08
color blue, dir08
show spheres, pEnd08
set sphere_scale, 0.25, pEnd08
color blue, pEnd08
label pEnd08, "8"
pseudoatom pEnd09, pos=[14.824,33.581,32.196]
distance dir09, p0, pEnd09
color blue, dir09
show spheres, pEnd09
set sphere_scale, 0.25, pEnd09
color blue, pEnd09
label pEnd09, "9"
pseudoatom pEnd10, pos=[18.894,21.165,34.215]
distance dir10, p0, pEnd10
color blue, dir10
show spheres, pEnd10
set sphere_scale, 0.25, pEnd10
color blue, pEnd10
label pEnd10, "10"
pseudoatom pEnd11, pos=[20.615,23.458,25.276]
distance dir11, p0, pEnd11
color blue, dir11
show spheres, pEnd11
set sphere_scale, 0.25, pEnd11
color blue, pEnd11
label pEnd11, "11"
pseudoatom pEnd12, pos=[18.974,24.050,27.632]
distance dir12, p0, pEnd12
color blue, dir12
show spheres, pEnd12
set sphere_scale, 0.25, pEnd12
color blue, pEnd12
label pEnd12, "12"
pseudoatom pEnd13, pos=[18.788,22.118,31.396]
distance dir13, p0, pEnd13
color blue, dir13
show spheres, pEnd13
set sphere_scale, 0.25, pEnd13
color blue, pEnd13
label pEnd13, "13"
pseudoatom pEnd14, pos=[54.703,23.302,49.740]
distance dir14, p0, pEnd14
color red, dir14
show spheres, pEnd14
set sphere_scale, 0.25, pEnd14
color red, pEnd14
label pEnd14, "14"
pseudoatom pEnd15, pos=[16.979,24.389,34.140]
distance dir15, p0, pEnd15
color blue, dir15
show spheres, pEnd15
set sphere_scale, 0.25, pEnd15
color blue, pEnd15
label pEnd15, "15"
pseudoatom pEnd16, pos=[24.253,18.789,25.883]
distance dir16, p0, pEnd16
color blue, dir16
show spheres, pEnd16
set sphere_scale, 0.25, pEnd16
color blue, pEnd16
label pEnd16, "16"
pseudoatom pEnd17, pos=[18.433,21.879,34.013]
distance dir17, p0, pEnd17
color blue, dir17
show spheres, pEnd17
set sphere_scale, 0.25, pEnd17
color blue, pEnd17
label pEnd17, "17"
pseudoatom pEnd18, pos=[19.632,21.922,29.221]
distance dir18, p0, pEnd18
color blue, dir18
show spheres, pEnd18
set sphere_scale, 0.25, pEnd18
color blue, pEnd18
label pEnd18, "18"
pseudoatom pEnd19, pos=[19.127,24.491,26.812]
distance dir19, p0, pEnd19
color blue, dir19
show spheres, pEnd19
set sphere_scale, 0.25, pEnd19
color blue, pEnd19
label pEnd19, "19"
pseudoatom pEnd20, pos=[23.076,16.962,40.738]
distance dir20, p0, pEnd20
color blue, dir20
show spheres, pEnd20
set sphere_scale, 0.25, pEnd20
color blue, pEnd20
label pEnd20, "20"
set label_distance_digits, 0
hide labels, dir*
