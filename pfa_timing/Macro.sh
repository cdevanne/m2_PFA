rm datas.txt

root -q -b -l rootFiles/300events5cm30pi-10kaon0.root 'Macro.C(5.0)'
root -q -b -l rootFiles/300events10cm30pi-10kaon0.root  'Macro.C(10.0)'
root -q -b -l rootFiles/300events15cm30pi-10kaon0.root  'Macro.C(15.0)'
root -q -b -l rootFiles/300events20cm30pi-10kaon0.root  'Macro.C(20.0)'
root -q -b -l rootFiles/300events25cm30pi-10kaon0.root  'Macro.C(25.0)'
root -q -b -l rootFiles/300events30cm30pi-10kaon0.root  'Macro.C(30.0)'

root  MacroAnalysis2.C




# rm datas.txt

# root -q -b -l rootFiles/5m100events5cm30pi-10kaon0.root 'Macro.C(5.0)'
# root -q -b -l rootFiles/5m100events10cm30pi-10kaon0.root  'Macro.C(10.0)'
# root -q -b -l rootFiles/5m100events15cm30pi-10kaon0.root  'Macro.C(15.0)'
# root -q -b -l rootFiles/5m100events20cm30pi-10kaon0.root  'Macro.C(20.0)'
# root -q -b -l rootFiles/5m100events25cm30pi-10kaon0.root  'Macro.C(25.0)'
# root -q -b -l rootFiles/5m100events30cm30pi-10kaon0.root  'Macro.C(30.0)'

# root  MacroAnalysis.C