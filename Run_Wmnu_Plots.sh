
ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt","RunF",80,0,200,"met_pt (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("mu_eta[0]","RunF",50,-2.5,2.5,"mu_eta")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("mu_pt[0]","RunF",80,0,200,"mu_pt (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("el_eta[0]","RunF",50,-2.5,2.5,"el_eta")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("el_pt[0]","RunF",80,0,200,"el_pt (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_phi","RunF",50,0,3.14,"met_phi")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("mu_relIso[0]","RunF",100,0,0.5,"mu_relIso")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("el_relIso[0]","RunF",100,0,0.5,"mu_relIso")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("npv","RunF",51,-0.5,50.5,"npv")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("njets","RunF",6,0,6,"njets")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt","RunF",80,0,200,"puppi_pt (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_phi","RunF",50,0,3.14,"puppi_phi")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_smeared","RunF",80,0,200,"met_pt_smeared (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT","RunF",80,0,200,"MT (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi","RunF",80,0,200,"MTpuppi (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_smeared","RunF",80,0,200,"MT_smeared (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut","RunF",80,0,800,"Ut (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr","RunF",80,-250,250,"Utr (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol","RunF",80,-250,100,"Ucol (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi","RunF",80,0,800,"Ut_puppi (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi","RunF",80,-250,250,"Utr_puppi (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi","RunF",80,-250,100,"Ucol_puppi (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared","RunF",80,0,800,"Ut_smeared (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared","RunF",80,-250,250,"Utr_smeared (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared","RunF",80,-250,100,"Ucol_smeared (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("PFmet_pt","RunF",80,0,200,"PFmet_pt (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii

ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_PF","RunF",80,0,200,"MT_PF (GeV)")'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!stop!!!!!!!!!!!!!!!!!!!!!!!!1
echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!stop!!!!!!!!!!!!!!!!!!!!!!!!1

echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!stop!!!!!!!!!!!!!!!!!!!!!!!!1

echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!stop!!!!!!!!!!!!!!!!!!!!!!!!1

### syst unc

###normal met
ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetEnUp","RunF",80,0,200,"met_pt_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetEnDown","RunF",80,0,200,"met_pt_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_UnclusteredEnUp","RunF",80,0,200,"met_pt_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_UnclusteredEnDown","RunF",80,0,200,"met_pt_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetResUp","RunF",80,0,200,"met_pt_JetResUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetResDown","RunF",80,0,200,"met_pt_JetResDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



###puppi met


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_JetEnUp","RunF",80,0,200,"puppi_pt_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_JetEnDown","RunF",80,0,200,"puppi_pt_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_UnclusteredEnUp","RunF",80,0,200,"puppi_pt_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_UnclusteredEnDown","RunF",80,0,200,"puppi_pt_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_JetResUp","RunF",80,0,200,"puppi_pt_JetResUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("puppi_pt_JetResDown","RunF",80,0,200,"puppi_pt_JetResDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




###smeared met
ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetEnUp_smeared","RunF",80,0,200,"met_pt_JetEnUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetEnDown_smeared","RunF",80,0,200,"met_pt_JetEnDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_UnclusteredEnUp_smeared","RunF",80,0,200,"met_pt_UnclusteredEnUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_UnclusteredEnDown_smeared","RunF",80,0,200,"met_pt_UnclusteredEnDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetResUp_smeared","RunF",80,0,200,"met_pt_JetResUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("met_pt_JetResDown_smeared","RunF",80,0,200,"met_pt_JetResDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


###syst for MT

###normal met
ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetEnUp","RunF",80,0,200,"MT_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetEnDown","RunF",80,0,200,"MT_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_UnclusteredEnUp","RunF",80,0,200,"MT_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_UnclusteredEnDown","RunF",80,0,200,"MT_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetResUp","RunF",80,0,200,"MT_JetResUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetResDown","RunF",80,0,200,"MT_JetResDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



###puppi met


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_JetEnUp","RunF",80,0,200,"MTpuppi_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_JetEnDown","RunF",80,0,200,"MTpuppi_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_UnclusteredEnUp","RunF",80,0,200,"MTpuppi_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_UnclusteredEnDown","RunF",80,0,200,"MTpuppi_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_JetResUp","RunF",80,0,200,"MTpuppi_JetResUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MTpuppi_JetResDown","RunF",80,0,200,"MTpuppi_JetResDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


###smeared met
ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetEnUp_smeared","RunF",80,0,200,"MT_JetEnUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetEnDown_smeared","RunF",80,0,200,"MT_JetEnDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_UnclusteredEnUp_smeared","RunF",80,0,200,"MT_UnclusteredEnUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_UnclusteredEnDown_smeared","RunF",80,0,200,"MT_UnclusteredEnDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetResUp_smeared","RunF",80,0,200,"MT_JetResUp_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_JetResDown_smeared","RunF",80,0,200,"MT_JetResDown_smeared",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


#################### Ut systematick

ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_JetResDown","RunF",80,0,800,"Ut_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_JetResDown","RunF",80,-250,250,"Utr_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_JetResDown","RunF",80,-250,100,"Ucol_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_JetResDown","RunF",80,0,800,"Ut_puppi_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_JetResDown","RunF",80,-250,250,"Utr_puppi_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_JetResDown","RunF",80,-250,100,"Ucol_puppi_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_JetResUp","RunF",80,0,800,"Ut_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_JetResUp","RunF",80,-250,250,"Utr_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_JetResUp","RunF",80,-250,100,"Ucol_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_JetResUp","RunF",80,0,800,"Ut_puppi_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_JetResUp","RunF",80,-250,250,"Utr_puppi_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_JetResUp","RunF",80,-250,100,"Ucol_puppi_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_JetEnUp","RunF",80,0,800,"Ut_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_JetEnUp","RunF",80,-250,250,"Utr_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_JetEnUp","RunF",80,-250,100,"Ucol_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_JetEnUp","RunF",80,0,800,"Ut_puppi_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_JetEnUp","RunF",80,-250,250,"Utr_puppi_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_JetEnUp","RunF",80,-250,100,"Ucol_puppi_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_JetEnDown","RunF",80,0,800,"Ut_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_JetEnDown","RunF",80,-250,250,"Utr_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_JetEnDown","RunF",80,-250,100,"Ucol_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_JetEnDown","RunF",80,0,800,"Ut_puppi_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_JetEnDown","RunF",80,-250,250,"Utr_puppi_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_JetEnDown","RunF",80,-250,100,"Ucol_puppi_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_UnclusteredEnDown","RunF",80,0,800,"Ut_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_UnclusteredEnDown","RunF",80,-250,250,"Utr_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_UnclusteredEnDown","RunF",80,-250,100,"Ucol_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_UnclusteredEnDown","RunF",80,0,800,"Ut_puppi_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_UnclusteredEnDown","RunF",80,-250,250,"Utr_puppi_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_UnclusteredEnDown","RunF",80,-250,100,"Ucol_puppi_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_UnclusteredEnUp","RunF",80,0,800,"Ut_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_UnclusteredEnUp","RunF",80,-250,250,"Utr_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_UnclusteredEnUp","RunF",80,-250,100,"Ucol_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_puppi_UnclusteredEnUp","RunF",80,0,800,"Ut_puppi_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_puppi_UnclusteredEnUp","RunF",80,-250,250,"Utr_puppi_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_puppi_UnclusteredEnUp","RunF",80,-250,100,"Ucol_puppi_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


#### smeared U



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_JetResDown","RunF",80,0,800,"Ut_smeared_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_JetResDown","RunF",80,-250,250,"Utr_smeared_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_JetResDown","RunF",80,-250,100,"Ucol_smeared_JetResDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii






ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_JetResUp","RunF",80,0,800,"Ut_smeared_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_JetResUp","RunF",80,-250,250,"Utr_smeared_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_JetResUp","RunF",80,-250,100,"Ucol_smeared_JetResUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii






ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_JetEnUp","RunF",80,0,800,"Ut_smeared_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_JetEnUp","RunF",80,-250,250,"Utr_smeared_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_JetEnUp","RunF",80,-250,100,"Ucol_smeared_JetEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_JetEnDown","RunF",80,0,800,"Ut_smeared_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_JetEnDown","RunF",80,-250,250,"Utr_smeared_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_JetEnDown","RunF",80,-250,100,"Ucol_smeared_JetEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii






ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_UnclusteredEnDown","RunF",80,0,800,"Ut_smeared_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_UnclusteredEnDown","RunF",80,-250,250,"Utr_smeared_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_UnclusteredEnDown","RunF",80,-250,100,"Ucol_smeared_UnclusteredEnDown (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ut_smeared_UnclusteredEnUp","RunF",80,0,800,"Ut_smeared_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Utr_smeared_UnclusteredEnUp","RunF",80,-250,250,"Utr_smeared_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("Ucol_smeared_UnclusteredEnUp","RunF",80,-250,100,"Ucol_smeared_UnclusteredEnUp (GeV)",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii




# PF met!!!!!



ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_PF_JetEnUp","RunF",80,0,200,"MT_PF_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_PF_JetEnDown","RunF",80,0,200,"MT_PF_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_PF_UnclusteredEnUp","RunF",80,0,200,"MT_PF_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("MT_PF_UnclusteredEnDown","RunF",80,0,200,"MT_PF_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii 

ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("PFmet_pt_JetEnUp","RunF",80,0,200,"PFmet_pt_JetEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("PFmet_pt_JetEnDown","RunF",80,0,200,"PFmet_pt_JetEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("PFmet_pt_UnclusteredEnUp","RunF",80,0,200,"PFmet_pt_UnclusteredEnUp",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii


ii=$((1 + RANDOM % 10000))
cat bss2 > Jobs/jobWmnu${ii}.sh
echo root -l -b -q "'"'PlotMy2.C("PFmet_pt_UnclusteredEnDown","RunF",80,0,200,"PFmet_pt_UnclusteredEnDown",true)'"'" >> Jobs/jobWmnu${ii}.sh
chmod u+x Jobs/jobWmnu${ii}.sh
./HTC_submit.sh  Jobs/jobWmnu${ii}.sh $ii



