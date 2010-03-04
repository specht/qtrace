e = 0.0
system("mkdir out")
while e <= 1.0 do
    command = "../../qtrace --label \"15N (#{sprintf('%1.4f', e)})\" --maxCharge 2 --xhtmlOutputTarget ./out/out-#{sprintf('%1.4f', e)}.svg --spectraFiles ../../scan-1.mzml --peptides LYPGGSFDPLGLADDPDTFAELK"
    system(command)
    e += 0.01
end
