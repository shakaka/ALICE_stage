

export PRINTER=irc-phenrv

export WORKDIRAFTER=/projet/AFTER/Cynthia

# User specific aliases and functions
alias l='ls -ltr'
alias win='rdesktop -k en-us -g 90% ipnnts4'
alias winbdd1='rdesktop -k en-us -a 16 -g 90% ipnbdd1.ipno.in2p3.fr'
alias xterm='xterm -sb -title ipno -name ipno' 

# ssh 
alias hermes='ssh -X -l cynthia geordi.desy.de' 
alias frascati='ssh -t -l cynthia -X axcalc.lnf.infn.it ssh -X ed22lf' 
alias windows='rdesktop -k en-us -r sound -u hadjida -a 16 -g 98% nants1 &'
alias easyfind='find ./ -name "\!:1" -exec grep -i \!:2 {} \; -print'
alias ccin2p3new='ssh -X -l hadjida ccali.in2p3.fr' 
alias ccin2p3='ssh -X -l hadjida ccali.in2p3.fr'
alias cern='ssh -X -l cynthia lxplus.cern.ch'
alias subatech='ssh -X -l hadjida nangate.in2p3.fr'
alias ipno='ssh -X -Y ipninter.in2p3.fr'
alias grid='ssh -X ipngrid01.in2p3.fr'

#aliroot
alias alirootv='source $WORKDIR/code/local/script/setup-alice-ubuntu.sh'
alias alirootvgrid='source $WORKDIR/code/local/script/SetEnvIPNGrid.sh '

# grille de calcul -GRIF-
## create a proxy
alias grid_proxyinit="voms-proxy-init -voms vo.ipno.in2p3.fr --valid 24:00" 
alias grid_myproxyinit="myproxy-init -s myproxy.grif.fr -d -n -t 24 -c 800"
alias grid_proxyinfo="voms-proxy-info -all"
alias grid_myproxyinfo="myproxy-info -s myproxy.grif.fr -d"
alias grid_proxyend="voms-proxy-destroy"
alias grid_myproxyend="myproxy-destroy"
## submit a job and get status, output, info
# wms ipngrid28 faster but do not work yet with InputData
#alias grid_jobsubmit="glite-wms-job-submit -a -e https://ipngrid28.in2p3.fr:7443/glite_wms_wmproxy_server -o"
#alias grid_joblistmatch="glite-wms-job-list-match -a -e https://ipngrid28.in2p3.fr:7443/glite_wms_wmproxy_server"
#alias grid_jobsubmit="glite-wms-job-submit -a -e https://node27.datagrid.cea.fr:7443/glite_wms_wmproxy_server -o"
#alias grid_joblistmatch="glite-wms-job-list-match -a -e https://node27.datagrid.cea.fr:7443/glite_wms_wmproxy_server"
alias grid_jobsubmit="glite-wms-job-submit  -a  -e https://marwms.in2p3.fr:7443/glite_wms_wmproxy_server -o"
alias grid_joblistmatch="glite-wms-job-list-match -a -e https://marwms.in2p3.fr:7443/glite_wms_wmproxy_server"
alias grid_jobstatus="glite-wms-job-status -v 0 -i"
alias grid_jobcancel="glite-wms-job-cancel -i" 
alias grid_joboutput="glite-wms-job-output -i"
alias grid_jobinfo="glite-wms-job-info --jdl -i"
## get ipno grid info on CE and SE
alias grid_lcginfo="lcg-infosites --vo vo.ipno.in2p3.fr all"

## Use SE
export LFC_HOME=/grid/vo.ipno.in2p3.fr/cynthia
#LCG command
alias se_ls="lcg-ls "
alias se_cp="lcg-cp "
alias se_cr="lcg-cr "
alias se_lr="lcg-lr "
alias se_rep="lcg-rep "
alias se_del="lcg-del -a "
#to be used to remove recursively a directory 
#se_rec_del $LFC_HOME/$file_path
alias se_rec_del="lcg-rec-del -vo vo.ipno.in2p3.fr -cp "
#guid
alias se_lg="lcg-lg "
#LFC command
alias se_mkdir="lfc-mkdir "
