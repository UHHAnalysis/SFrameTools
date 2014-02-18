#!/usr/bin/env python

import datablock
import math
from optparse import OptionParser
import os
import re
import ROOT
import shutil
import sys
import time
import xmlparser

parser = OptionParser()
parser.add_option("-c", "--cfg", dest="configxml", help="Input XML Config File")
parser.add_option("-j", "--jobname", dest="jobname", default="", help="Job Name")
parser.add_option("-k", "--kill", dest="kill",  default="none", help="Kill Jobs: all, 1-5, or 1,3,5,9")
parser.add_option("-n", "--numjobs", dest="numjobs", type="int", default=0, help="Number of Jobs")
parser.add_option("-s", "--submit", dest="submit",  default="none", help="Submit Jobs: all, 1-5, or 1,3,5,9")

parser.add_option("--clobber", action="store_true", dest="clobber", default=False, help="Overwrite Job Directory")
parser.add_option("--create", action="store_true", dest="create", default=False, help="Create job and configuration files.")
parser.add_option("--status", action="store_true", dest="status", default=False, help="Get job status")
parser.add_option("--retar", action="store_true", dest="retar", default=False, help="Recreate the job tarball.")
(options, args) = parser.parse_args()
if options.jobname == "": options.jobname=options.configxml.strip(".xml")
if options.jobname == "" and options.configxml == "":
    print "ERROR: Please provide either a configuration file or job directory"
    exit(1)

def makejoblist(joblist):
    newjoblist=[]
    joblistarray=joblist.split(",")
    for job in joblistarray:
        if job.find("-")!=-1:
            for i in range(int(job.split('-')[0]),int(job.split('-')[1])+1):
                newjoblist.append(i)
        else: newjoblist.append(int(job))
    return newjoblist

def makestringlist(joblist):
    stringlist=""
    for job in joblist: stringlist += str(job)+","
    stringlist=stringlist[:-1]
    return stringlist

def getoutputfilenames(configfile):
    rawxmlfile = open(configfile).read()
    typelist = xmlparser.parse(rawxmlfile,"Type")
    versionlist = xmlparser.parse(rawxmlfile,"Version")
    namelist = xmlparser.parse(rawxmlfile,"Name")
    rootfilenamelist = ""
    for name in namelist:
        if name.find("Cycle") != -1: cyclename=name
    for type,version in zip(typelist,versionlist):
        rootfilenamelist += cyclename+"."+type+"."+version+".root, "
    return rootfilenamelist[:-2]

def createcondortxt(jobname, jobnumber):
    rootfiles = getoutputfilenames(jobname+"/xml/"+jobname+"_"+str(jobnumber)+".xml")
    os.chdir(jobname+"/configs")
    condorfile = open(jobname+"_"+str(jobnumber)+".txt", 'w')
    print >> condorfile, """universe = vanilla
Executable = %s/configs/%s_%d.sh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
InitialDir = %s/results
Transfer_Input_Files = ../configs/%s.tgz
Transfer_Output_Files = %s
Output = ../logs/%s_%d.stdout
Error = ../logs/%s_%d.stderr
Log = ../logs/%s_%d.log
notify_user = ${LOGNAME}@FNAL.GOV
Arguments = %s 0
Queue 1""" %(jobname, jobname, jobnumber, jobname, jobname, rootfiles.replace(".root","."+str(jobnumber)+".root"), jobname, jobnumber, jobname, jobnumber, jobname, jobnumber, os.getcwd())
    condorfile.close()
    os.chdir("../..")

def createcondorscript(jobname, jobnumber,eosstatusdir):
    os.chdir(jobname+"/configs")
    scriptname = jobname+"_"+str(jobnumber)+".sh"
    scriptfile = open(scriptname, 'w')
    print >> scriptfile, """#!/bin/bash
WORKINGDIR=$PWD
STATUSFILE=%s/%s_%d.status
echo 'Configuring' >& $STATUSFILE
TARNAME=`/bin/ls *.tgz`
tar -xzf $TARNAME
SFRAMEDIR=`find . -type f -name fullsetup.sh | xargs dirname`
cd $SFRAMEDIR
eval `scramv1 runtime -sh`
source fullsetup.sh
cd $WORKINGDIR
ANALYSISDIR=CMSSW`echo ${1##*CMSSW}`/../../
cd $ANALYSISDIR
cp %s/xml/%s_%d.xml .
echo 'Running' >& $STATUSFILE
sframe_main %s_%d.xml
for filename in `/bin/ls *.root`; do
    newfilename=`echo $filename | sed 's|.root|.%d.root|'`
    mv $filename $newfilename
done
mv *.root $WORKINGDIR
echo 'Done' >& $STATUSFILE""" %(eosstatusdir,jobname, jobnumber, jobname, jobname, jobnumber, jobname, jobnumber, jobnumber)
    os.chmod(scriptname, 493) #493==755 in python chmod
    os.chdir("../..")

def begindatablock(mydatablock):
    return '<InputData Lumi="'+mydatablock.blocklumi+'" NEventsMax="'+mydatablock.neventsmax+'" Type="'+mydatablock.type+'" Version="'+mydatablock.version+'" Cacheable="'+mydatablock.cacheable+'">\n'

def enddatablack(mydatablock,indent):
    returnstring = indent+"  "+'<InputTree Name="AnalysisTree" />\n'
    if len(mydatablock.namelist)>1: returnstring += indent+"  "+'<OutputTree Name="AnalysisTree" />\n'
    returnstring += indent+"</InputData>\n"
    return returnstring

def makedatablocks(xmlfile):
    input=xmlfile
    datablocklist=[]
    while input.find("<InputData ") != -1:
        datablock = input[input.find("<InputData "):input.find("</InputData>")+12]
        datablocklist.append(datablock)
        input = input[input.find("</InputData>")+12:]
    return datablocklist

def resolveentities(input):
    begincomment = input.find("<!--")
    endcomment = input.find("-->")
    while begincomment!=-1 and endcomment!=-1:
        input = input[:begincomment]+input[endcomment+3:]
        begincomment = input.find("<!--")
        endcomment = input.find("-->")
    inputlist = input.split("\n")
    entitydict = {}
    for line in inputlist:
        if line.count("ENTITY")>0:
            entityname = line[line.find("ENTITY ")+7:line.find(" SYSTEM")]
            entityfile = line[line.find('"')+1:line.rfind('"')]
            entityfiles = open(entityfile).read()
            entitytemp = {entityname: entityfiles}
            entitydict.update(entitytemp)
    for entity in entitydict: input=input.replace("&"+entity+";", entitydict[entity])
    xmldump=open("dump.xml",'w');xmldump.write(input);xmldump.close()
    return input

def createxmlfiles(configfile, jobname, numjobs):
    rawxmlfile = open(configfile).read()
    xmlfile = resolveentities(rawxmlfile)
    xmldatablocks = makedatablocks(xmlfile)
    totalfilelist = xmlparser.parse(xmlfile,"FileName")
    datablocklist = []
    for xmldatablock in xmldatablocks:
        type = xmlparser.parse(xmldatablock,"Type")[0]
        version = xmlparser.parse(xmldatablock,"Version")[0]
        maxevents = xmlparser.parse(xmldatablock,"NEventsMax")[0]
        cacheable = xmlparser.parse(xmldatablock,"Cacheable")[0]
        filelist = xmlparser.parse(xmldatablock,"FileName")
        lumi = xmlparser.parse(xmldatablock,"Lumi")
        blocklumi = lumi.pop(0)
        namelist = xmlparser.parse(xmldatablock,"Name")
        mydatablock=datablock.datablock(blocklumi, filelist, lumi, namelist, type, version, maxevents, cacheable)
        datablocklist.append(mydatablock)
    if numjobs==0:
        totalfilelist = xmlparser.parse(xmlfile,"FileName")
        numjobs=len(totalfilelist)
        options.numjobs=len(totalfilelist)
    numjoblist=[]
    for i in range(options.numjobs): numjoblist.append(0)
    jobindex = 0
    for i in range(len(totalfilelist)):
        if jobindex>=options.numjobs: jobindex-=options.numjobs
        numjoblist[jobindex]+=1
        jobindex+=1
    datablocknumber = 0
    blockindex = 0
    for jobnumber in range(options.numjobs):
        blockindex,datablocknumber = createxmlfile(xmlfile, jobname, jobnumber+1, datablocklist, datablocknumber, blockindex, numjoblist[jobnumber])

def createxmlfile(infile, jobname, jobnumber, datablocklist, datablocknumber, blockindex, numfiles):
    filename = jobname+"_"+str(jobnumber)+".xml"
    os.chdir(jobname+"/xml")
    frontend = infile[:infile.find("<InputData ")]
    indent = frontend[frontend.rfind("\n")+1:]
    backend = infile[infile.rfind("</InputData>")+13:]

    if len(datablocklist)>0:
        inputfilestring = begindatablock(datablocklist[datablocknumber])
        filelist=datablocklist[datablocknumber].filelist
        lumilist=datablocklist[datablocknumber].lumilist
    else:
        print "Error: No DataBlocks Found!\n"
        exit(2)

    while numfiles>0:
        if blockindex == len(filelist):
            inputfilestring += enddatablack(datablocklist[datablocknumber],indent)
            datablocknumber += 1
            if datablocknumber < len(datablocklist):
                inputfilestring += indent+begindatablock(datablocklist[datablocknumber])
                filelist=datablocklist[datablocknumber].filelist
                lumilist=datablocklist[datablocknumber].lumilist
                blockindex=0
        else:
            inputfilestring += indent+"    "+'<In FileName="'+filelist[blockindex]+'" Lumi="'+lumilist[blockindex]+'"/>\n'
            blockindex += 1
            numfiles -= 1
    inputfilestring += enddatablack(datablocklist[datablocknumber],indent)

    if blockindex == len(filelist):
        datablocknumber += 1
        blockindex = 0

    finalxmlfile = frontend + inputfilestring + backend
    outputdirectory = re.search('OutputDirectory="[^ ]*"',finalxmlfile)
    if not outputdirectory: print "Error: OutputDirectory not found in "+infile
    finalxmlfile = finalxmlfile.replace(outputdirectory.group(0),'OutputDirectory="./"')
    outfile = open(filename,'w')
    outfile.write(finalxmlfile)
    outfile.close()
    os.chdir("../..")
    return blockindex,datablocknumber

def getcondornode(jobid):
    condorstatus=os.popen("condor_q -submitter $USER").readlines()
    for linenumber in range(len(condorstatus)):
        if condorstatus[linenumber].find(jobid) != -1:
            for nodeline in reversed(range(0,linenumber)):
                if condorstatus[nodeline].find("Submitter:") != -1:
                    statuslist = condorstatus[nodeline].split(" ")
                    nodeindex = len(statuslist)-1
                    return statuslist[nodeindex].strip("\n")
    print "Error: Condor node with jobid: "+jobid+" not found!"
    exit(3)

def checklog(jobname, jobnumber):
    errorinfo = os.popen("egrep -i 'exit|break|exceed|error|traceback|aborted' "+jobname+"/logs/"+jobname+"_"+str(jobnumber)+".log | tr '\n' ', '").readline().strip("\n")
    errorline = errorinfo[errorinfo.find(")")+2:]
    return errorline

def checkstdout(jobname, jobnumber):
    errorinfo = os.popen('egrep -i "exit|break|exceed|error|traceback|aborted|E R R O R|find tree AnalysisTree|fatal" '+jobname+"/logs/"+jobname+"_"+str(jobnumber)+".stdout | tr '\n' ', '").readline().strip("\n")
    if errorinfo.find(":") != -1: errorline = errorinfo.split(":")[1]
    else: errorline=errorinfo
    return errorline

def getjobinfo(jobname,jobnumber,resubmitjobs):
    outputfiles = os.popen("grep Transfer_Output_Files "+jobname+"/configs/"+jobname+"_"+str(jobnumber)+".txt").readline().strip("\n")
    outputfiles = outputfiles.split(" ")[2:]
    jobinfo=""
    for file in outputfiles:
        file=file.strip(",")
        filepath = jobname+"/results/"+file
        jobstatus = open(options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status").read().strip("\n")
        if os.path.isfile(filepath):
            if os.path.getsize(filepath)==0:
                if jobstatus=="Done": jobinfo += " Root File "+file+" Empty!"
            else:
                rootfile = ROOT.TFile.Open(filepath)
                if rootfile==None:
                    resubmitjobs.append(jobnumber)
                    jobinfo += " Root File "+file+" is a Null Pointer!"
                elif rootfile.IsZombie():
                    resubmitjobs.append(jobnumber)
                    jobinfo += " Root File "+file+" is Zombie!"
                elif not rootfile.IsZombie():
                    analysistree = rootfile.Get("AnalysisTree")
                    jobinfo += " Root File "+file+" is Valid: "+str(int(analysistree.GetEntries()))+" Events."
                else:
                    hist = rootfile.Get("nprocessed")
                    jobinfo += " Root File "+file+" is Valid: "+str(int(hist.GetEntries()))+" Events."
                    if file.find("PostSelection") == -1: jobinfo += " Warning No AnalysisTree Found in "+file
        else:
            jobinfo += " Output file "+file+" is not found!"
        if jobstatus=="Done":
            logerror = checklog(jobname, jobnumber)
            stdouterror = checkstdout(jobname, jobnumber)
            if (logerror != "" or stdouterror != "") and resubmitjobs.count(jobnumber)<1: resubmitjobs.append(jobnumber)
            if logerror != "": jobinfo += " "+logerror
            if stdouterror != "": jobinfo += " "+stdouterror
    return jobinfo

if not options.create and options.submit=="none" and options.kill=="none" and not options.status:
    print "ERROR: Must either create, submit jobs, kill, or check the status of jobs"

workingdir=os.getcwd()
cmsswbase=os.getenv("CMSSW_BASE")
username=os.getenv("USER")
currentnode=os.getenv("HOST")
eosstatusdir="/eos/uscms/store/user/"+username+"/BSFrameStatus/"+options.jobname+"/"

if options.create:
    if not os.path.isdir(options.jobname) and not options.clobber: os.mkdir(options.jobname)
    elif os.path.isdir(options.jobname) and not options.clobber:
        print "ERROR: Job directory "+options.jobname+" exists!\nPlease remove job directory or enable --clobber option.."
        exit(4)
    elif options.clobber:
        if os.path.isdir(options.jobname): shutil.rmtree(options.jobname)
        os.mkdir(options.jobname)
    os.mkdir(options.jobname+"/configs")
    os.mkdir(options.jobname+"/logs")
    os.mkdir(options.jobname+"/results")
    os.mkdir(options.jobname+"/status")
    os.mkdir(options.jobname+"/xml")

    if not os.path.isdir("/eos/uscms/store/user/"+username+"/BSFrameStatus/"): os.mkdir("/eos/uscms/store/user/"+username+"/BSFrameStatus/")
    if not os.path.isdir(eosstatusdir): os.mkdir(eosstatusdir)
    else:
        shutil.rmtree(eosstatusdir)
        os.mkdir(eosstatusdir)

    createxmlfiles(options.configxml,options.jobname,options.numjobs)
    for jobnumber in range(1,options.numjobs+1):
        createcondortxt(options.jobname,jobnumber)
        createcondorscript(options.jobname,jobnumber,eosstatusdir)
        os.system("echo 'Created' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

    os.chdir(cmsswbase+"/..")
    tarball = options.jobname+".tgz"
    target = os.popen("echo ${CMSSW_BASE##*/}").readline().strip("\n")+"/"
    print "Creating tarball of "+target+" area."
    os.system("tar -czf "+tarball+" "+target+" --exclude='*.root' --exclude='*.tgz'")
    os.system("mv "+tarball+" "+workingdir+"/"+options.jobname+"/configs")
    os.chdir(workingdir)

if not os.path.isdir(options.jobname):
    print "ERROR: Job directory "+options.jobname+" does not exist!\nPlease create job with bsframe.py -c myconfig.xml --create."
    exit(5)

if options.retar:
    if options.create: print "There is no point in creating a task and then recreating the tarball."
    os.chdir(cmsswbase+"/..")
    tarball = options.jobname+".tgz"
    target = os.popen("echo ${CMSSW_BASE##*/}").readline().strip("\n")+"/"
    print "Creating tarball of "+target+" area."
    os.system("tar -czf "+tarball+" "+target+" --exclude='*.root' --exclude='*.tgz'")
    os.system("mv "+tarball+" "+workingdir+"/"+options.jobname+"/configs")
    os.chdir(workingdir)

if options.kill!="none":
    joblist=[]
    if options.kill=="all":
        if options.numjobs==0: options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip('\n'))
        joblist=range(1,options.numjobs+1)
    else: joblist=makejoblist(options.kill)
    condorstatus=os.popen("condor_q -submitter $USER").read()
    for jobnumber in joblist:
        print "Killing job number: %d" %(jobnumber)
        logfile = os.popen("/bin/ls -rt "+options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log | tail -1").readline().strip('\n')
        jobid = os.popen("grep submitted "+logfile+" | tail -1 | awk '{print $2}'").readline().strip("\n()").split(".")[0]
        subnode = getcondornode(jobid)
        if condorstatus.find(jobid) == -1: print "Error! Job "+jobid+" not found."
        else:
            if subnode==currentnode: os.system("condor_rm "+jobid)
            else: os.system("rsh "+subnode+" condor_rm "+jobid)
            time.sleep(0.3)
            os.system("echo 'Killed' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

if options.submit!="none":
    joblist=[]
    if options.submit=="all":
        if options.numjobs==0: options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip('\n'))
        joblist=range(1,options.numjobs+1)
    else: joblist=makejoblist(options.submit)
    print "Submitting %d jobs" %(len(joblist))
    for jobnumber in joblist:
        print "Submitting job number: %d" %(jobnumber)
        subnum = int(os.popen("grep Arguments "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt | awk '{print $4}'").readline().strip('\n'))
        os.system("sed -i 's|Transfer_Output_Files = \(.*\)_"+str(subnum)+".root$|Transfer_Output_Files = \\1_"+str(subnum+1)+".root|' "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        os.system("sed -i 's|/configs "+str(subnum)+"$|/configs "+str(subnum+1)+"|' "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        if os.path.isfile(options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log"): os.system("/bin/rm "+options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log")
        os.system("condor_submit "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        os.system("echo 'Submitted' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

jobstatuslist=[]
resubmitjobs=[]
if (options.status):
    for statuslog in os.popen("/bin/ls "+eosstatusdir).readlines():
        statuslog=statuslog.strip("\n")
        filesize = int(os.popen("ls -l "+eosstatusdir+"/"+statuslog+" | awk '{print $5}'").readline().strip("\n"))
        eostimestamp = os.path.getmtime(eosstatusdir+"/"+statuslog)
        localtimestamp = os.path.getmtime(options.jobname+"/status/"+statuslog)
        if filesize>0 and eostimestamp>localtimestamp: os.system("/bin/cp "+eosstatusdir+"/"+statuslog+" "+options.jobname+"/status/")
    print "\nJob Status Summary for Task: ",options.jobname
    print "================================================================================"
    print "Job Number         Status             Additional Information"
    print "--------------------------------------------------------------------------------"
    whitespace="                                                                                "
    if options.numjobs==0: options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip("\n"))
    for jobnumber in range(1,options.numjobs+1):
        jobstatus=open(options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status").read().strip("\n")
        jobstatuslist.append(jobstatus)
        jobinfo=getjobinfo(options.jobname,jobnumber,resubmitjobs)
        print whitespace[:4]+str(jobnumber)+whitespace[:15-len(str(jobnumber))]+jobstatus+whitespace[:18-len(jobstatus)]+jobinfo
    print ""
    if jobstatuslist.count("Created")>0: print "There are "+str(jobstatuslist.count("Created"))+" Created Jobs"
    if jobstatuslist.count("Submitted")>0: print "There are "+str(jobstatuslist.count("Submitted"))+" Submitted Jobs"
    if jobstatuslist.count("Configuring")>0: print "There are "+str(jobstatuslist.count("Configuring"))+" Configuring Jobs"
    if jobstatuslist.count("Running")>0: print "There are "+str(jobstatuslist.count("Running"))+" Running Jobs"
    if jobstatuslist.count("Killed")>0: print "There are "+str(jobstatuslist.count("Killed"))+" Killed Jobs"
    if jobstatuslist.count("Done")>0: print "There are "+str(jobstatuslist.count("Done"))+" Done Jobs"
    if len(resubmitjobs)>0:
        print "\nThere are "+str(len(resubmitjobs))+" jobs with problems!!!!"
        print "Resubmit jobs: "+makestringlist(resubmitjobs)

print "\nBSFrame has completed execution!"
