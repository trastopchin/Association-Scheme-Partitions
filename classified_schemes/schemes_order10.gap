#
# computed by Izumi Miyamoto and Akihide Hanaki 
# using gap 3.4.4
#
as10 := 
[
# No. 1
[ [ 0,1,1,1,1,1,1,1,1,1 ],
  [ 1,0,1,1,1,1,1,1,1,1 ],
  [ 1,1,0,1,1,1,1,1,1,1 ],
  [ 1,1,1,0,1,1,1,1,1,1 ],
  [ 1,1,1,1,0,1,1,1,1,1 ],
  [ 1,1,1,1,1,0,1,1,1,1 ],
  [ 1,1,1,1,1,1,0,1,1,1 ],
  [ 1,1,1,1,1,1,1,0,1,1 ],
  [ 1,1,1,1,1,1,1,1,0,1 ],
  [ 1,1,1,1,1,1,1,1,1,0 ] ]
,
# No. 2
[ [ 0,1,2,2,2,2,2,2,2,2 ],
  [ 1,0,2,2,2,2,2,2,2,2 ],
  [ 2,2,0,1,2,2,2,2,2,2 ],
  [ 2,2,1,0,2,2,2,2,2,2 ],
  [ 2,2,2,2,0,1,2,2,2,2 ],
  [ 2,2,2,2,1,0,2,2,2,2 ],
  [ 2,2,2,2,2,2,0,1,2,2 ],
  [ 2,2,2,2,2,2,1,0,2,2 ],
  [ 2,2,2,2,2,2,2,2,0,1 ],
  [ 2,2,2,2,2,2,2,2,1,0 ] ]
,
# No. 3
[ [ 0,1,1,1,2,2,2,2,2,2 ],
  [ 1,0,2,2,1,1,2,2,2,2 ],
  [ 1,2,0,2,2,2,1,1,2,2 ],
  [ 1,2,2,0,2,2,2,2,1,1 ],
  [ 2,1,2,2,0,2,1,2,1,2 ],
  [ 2,1,2,2,2,0,2,1,2,1 ],
  [ 2,2,1,2,1,2,0,2,2,1 ],
  [ 2,2,1,2,2,1,2,0,1,2 ],
  [ 2,2,2,1,1,2,2,1,0,2 ],
  [ 2,2,2,1,2,1,1,2,2,0 ] ]
,
# No. 4
[ [ 0,1,1,1,1,2,2,2,2,2 ],
  [ 1,0,1,1,1,2,2,2,2,2 ],
  [ 1,1,0,1,1,2,2,2,2,2 ],
  [ 1,1,1,0,1,2,2,2,2,2 ],
  [ 1,1,1,1,0,2,2,2,2,2 ],
  [ 2,2,2,2,2,0,1,1,1,1 ],
  [ 2,2,2,2,2,1,0,1,1,1 ],
  [ 2,2,2,2,2,1,1,0,1,1 ],
  [ 2,2,2,2,2,1,1,1,0,1 ],
  [ 2,2,2,2,2,1,1,1,1,0 ] ]
,
# No. 5
[ [ 0,1,2,2,2,2,3,3,3,3 ],
  [ 1,0,2,2,2,2,3,3,3,3 ],
  [ 2,2,0,1,3,3,2,2,3,3 ],
  [ 2,2,1,0,3,3,2,2,3,3 ],
  [ 2,2,3,3,0,1,3,3,2,2 ],
  [ 2,2,3,3,1,0,3,3,2,2 ],
  [ 3,3,2,2,3,3,0,1,2,2 ],
  [ 3,3,2,2,3,3,1,0,2,2 ],
  [ 3,3,3,3,2,2,2,2,0,1 ],
  [ 3,3,3,3,2,2,2,2,1,0 ] ]
,
# No. 6
[ [ 0,1,2,2,2,2,3,3,3,3 ],
  [ 1,0,3,3,3,3,2,2,2,2 ],
  [ 2,3,0,2,2,2,1,3,3,3 ],
  [ 2,3,2,0,2,2,3,1,3,3 ],
  [ 2,3,2,2,0,2,3,3,1,3 ],
  [ 2,3,2,2,2,0,3,3,3,1 ],
  [ 3,2,1,3,3,3,0,2,2,2 ],
  [ 3,2,3,1,3,3,2,0,2,2 ],
  [ 3,2,3,3,1,3,2,2,0,2 ],
  [ 3,2,3,3,3,1,2,2,2,0 ] ]
,
# No. 7
[ [ 0,1,1,2,2,3,3,3,3,3 ],
  [ 1,0,2,1,2,3,3,3,3,3 ],
  [ 1,2,0,2,1,3,3,3,3,3 ],
  [ 2,1,2,0,1,3,3,3,3,3 ],
  [ 2,2,1,1,0,3,3,3,3,3 ],
  [ 3,3,3,3,3,0,1,1,2,2 ],
  [ 3,3,3,3,3,1,0,2,1,2 ],
  [ 3,3,3,3,3,1,2,0,2,1 ],
  [ 3,3,3,3,3,2,1,2,0,1 ],
  [ 3,3,3,3,3,2,2,1,1,0 ] ]
,
# No. 8
[ [ 0,1,2,3,4,5,5,5,5,5 ],
  [ 2,0,3,4,1,5,5,5,5,5 ],
  [ 1,4,0,2,3,5,5,5,5,5 ],
  [ 4,3,1,0,2,5,5,5,5,5 ],
  [ 3,2,4,1,0,5,5,5,5,5 ],
  [ 5,5,5,5,5,0,1,2,3,4 ],
  [ 5,5,5,5,5,2,0,3,4,1 ],
  [ 5,5,5,5,5,1,4,0,2,3 ],
  [ 5,5,5,5,5,4,3,1,0,2 ],
  [ 5,5,5,5,5,3,2,4,1,0 ] ]
,
# No. 9
[ [ 0,1,2,2,3,3,4,4,5,5 ],
  [ 1,0,3,3,2,2,5,5,4,4 ],
  [ 2,3,0,4,1,5,2,4,3,5 ],
  [ 2,3,4,0,5,1,4,2,5,3 ],
  [ 3,2,1,5,0,4,3,5,2,4 ],
  [ 3,2,5,1,4,0,5,3,4,2 ],
  [ 4,5,2,4,3,5,0,2,1,3 ],
  [ 4,5,4,2,5,3,2,0,3,1 ],
  [ 5,4,3,5,2,4,1,3,0,2 ],
  [ 5,4,5,3,4,2,3,1,2,0 ] ]
,
# No. 10
[ [ 0,1,2,2,3,3,4,4,5,5 ],
  [ 1,0,4,4,5,5,2,2,3,3 ],
  [ 2,5,0,3,2,3,4,5,1,4 ],
  [ 2,5,3,0,3,2,5,4,4,1 ],
  [ 3,4,2,3,0,2,5,1,5,4 ],
  [ 3,4,3,2,2,0,1,5,4,5 ],
  [ 5,2,5,4,4,1,0,3,2,3 ],
  [ 5,2,4,5,1,4,3,0,3,2 ],
  [ 4,3,1,5,4,5,2,3,0,2 ],
  [ 4,3,5,1,5,4,3,2,2,0 ] ]
,
# No. 11
[ [ 0,1,2,2,3,3,4,4,5,5 ],
  [ 1,0,2,2,3,3,4,4,5,5 ],
  [ 3,3,0,1,4,4,5,5,2,2 ],
  [ 3,3,1,0,4,4,5,5,2,2 ],
  [ 2,2,5,5,0,1,3,3,4,4 ],
  [ 2,2,5,5,1,0,3,3,4,4 ],
  [ 5,5,4,4,2,2,0,1,3,3 ],
  [ 5,5,4,4,2,2,1,0,3,3 ],
  [ 4,4,3,3,5,5,2,2,0,1 ],
  [ 4,4,3,3,5,5,2,2,1,0 ] ]
];
### end of file ###
