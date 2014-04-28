#!/usr/bin/env python

import datablock
import math
from optparse import OptionParser
import os
import re
import shutil
import sys
import time
import xmlparser

parser = OptionParser()
parser.add_option("-c", "--cfg", dest="configxml", default="", help="Input XML Config File")
parser.add_option("-j", "--jobname", dest="jobname", default="", help="Job Name")
parser.add_option("-k", "--kill", dest="kill",  default="", help="Kill Jobs: all, 1-5, or 1,3,5,9")
parser.add_option("-n", "--numjobs", dest="numjobs", type="int", default=0, help="Number of Jobs")
parser.add_option("-o", "--output", dest="output", default="", help="Output directory: The default is jobname/results")
parser.add_option("-s", "--submit", dest="submit",  default="", help="Submit Jobs: all, 1-5, or 1,3,5,9")
parser.add_option("-l", "--clean", dest="clean",  default="", help="Clean status, log, and result files: all, 1-5, 1,3,5,9")

parser.add_option("--clobber", action="store_true", dest="clobber", default=False, help="Overwrite Job Directory")
parser.add_option("--create", action="store_true", dest="create", default=False, help="Create job and configuration files.")
parser.add_option("--status", action="store_true", dest="status", default=False, help="Get job status")
parser.add_option("--retar", action="store_true", dest="retar", default=False, help="Recreate the job tarball.")
parser.add_option("--ttbargencut", action="store_true", dest="ttbargencut", default=False, help="Apply ttbar generator cut")
parser.add_option("--notar", action="store_true", dest="notar", default=False, help="Do not create tarball. For debugging configs.")

parser.add_option("--append", dest="append", default="", help="Append string to job name")
parser.add_option("--flavor", dest="flavor", default="", help="Apply flavor selection: bflavor, cflavor, lflavor")
parser.add_option("--pileupfile", dest="pileupfile", default="", help="Specify pileup file")
parser.add_option("--bjets", dest="bjets", default="", help="Apply bjet Systematic: up-bjets, down-bjets, up-ljets, down-ljets")
parser.add_option("--toptag", dest="tjets", default="", help="Apply toptag scale Systematic: up-mistag, down-mistag, up-toptag, down-toptag")
parser.add_option("--JEC", dest="jec", default="", help="Apply JEC Systematic: up or down")
parser.add_option("--JER", dest="jer", default="", help="Apply JER Systematic: up or down")
parser.add_option("--PDF", dest="pdf", default="", help="Apply PDF Systematics: CT10 or cteq66")
parser.add_option("--PDFDir", dest="pdfdir", default="", help="Location of PDF systematic files.")
parser.add_option("--filter", dest="filter", default="", help="Run only samples that pass filter.")
parser.add_option("--veto", dest="veto", default="", help="Remove samples that pass filter.")
(options, args) = parser.parse_args()

if options.jobname == "":
    options.jobname = options.configxml.strip(".xml")
    if options.ttbargencut: options.jobname += "_TTBar"
    if options.flavor != "": options.jobname += "_"+options.flavor
    if options.bjets != "": options.jobname += "_"+options.bjets
    if options.jec != "": options.jobname += "_JEC"+options.jec
    if options.jer != "": options.jobname += "_JER"+options.jer
    if options.pdf != "": options.jobname += "_"+options.pdf
    if options.append != "": options.jobname += "_"+options.append
if options.jobname == "" and options.configxml == "":
    print "ERROR: Please provide either a configuration file or job directory"
    exit(1)
if not os.path.isfile(options.configxml) and not os.path.isdir(options.jobname):
    if options.configxml!="": print "ERROR: "+options.configxml+" is not a valid file!"
    else: print "ERROR: "+options.jobname+" is not a valid directory!"
    exit(2)

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

def applypostfix(infile, append):
    postfix = re.search('PostFix="[^ ]+"',infile)
    if postfix: infile = infile.replace(postfix.group(0),postfix.group(0)[:-1]+"_"+str(append)+'"')
    elif infile.find('PostFix=""') != -1: infile = infile.replace('PostFix=""','PostFix="_'+str(append)+'"')
    else:
        cyclepos = infile[infile.find("<Cycle "):].find(">")+infile.find("<Cycle ")
        infile = infile[:cyclepos]+'PostFix="_'+str(append)+'" '+infile[cyclepos:]
    return infile

def additem(infile, name, value):
    location = re.search('Name="'+name+'" Value="[^ ]+"',infile)
    if location: infile = infile.replace(location.group(0),'Name="'+name+'" Value="'+str(value)+'"')
    else:
        configpos = infile.find("</UserConfig>")
        infile = infile[:configpos]+'  <Item Name="'+name+'" Value="'+str(value)+'" />\n    '+infile[configpos:]
    return infile

def applyttbargencut(infile):
    infile=applypostfix(infile,"0to700")
    infile=additem(infile,"ApplyMttbarGenCut","True")
    return infile

def applybjetsystematic(infile,bjets):
    infile=applypostfix(infile,bjets)
    infile=additem(infile,"BTaggingScaleFactors",bjets)
    return infile

def applyflavorselection(infile,flavor):
    infile=applypostfix(infile,flavor)
    infile=additem(infile,"ApplyFlavorSelection",flavor)
    return infile

def applyjesystematic(infile,jectype,jecdirection):
    infile=applypostfix(infile,jectype+jecdirection)
    infile=additem(infile,"SystematicUncertainty",jectype)
    infile=additem(infile,"SystematicVariation",jecdirection)
    return infile

def applypdfsystematics(infile, options, pdfindex):
    infile=applypostfix(infile, options.pdf+"_"+str(pdfindex))
    infile=additem(infile,"SystematicUncertainty","PDF")
    infile=additem(infile,"SystematicVariation","up")
    infile=additem(infile,"PDFName",options.pdf)
    if options.pdfdir != "": infile=additem(infile,"PDFWeightFilesDirectory",options.pdfdir)
    infile=additem(infile,"PDFIndex",pdfindex)
    return infile

def applytjetsystematic(infile,tjets):
    infile=applypostfix(infile,tjets)
    infile=additem(infile,"TopTaggingScaleFactors",tjets)
    return infile

def changepileupfile(infile,pileupfile):
    infile = additem(infile,"PU_Filename_Data",pileupfile)
    return infile

def getoutputfilenames(configfile):
    rawxmlfile = open(configfile).read()
    typelist = xmlparser.parse(rawxmlfile,"Type")
    versionlist = xmlparser.parse(rawxmlfile,"Version")
    namelist = xmlparser.parse(rawxmlfile,"Name")
    postfix = xmlparser.parse(rawxmlfile,"PostFix")
    rootfilenamelist = ""
    for name in namelist:
        if name.find("SelectionCycle") != -1: cyclename=name
    for type,version in zip(typelist,versionlist):
        rootfilenamelist += cyclename+"."+type+"."+version+postfix[0]+".root, "
    return rootfilenamelist[:-2]

def getinputfilenames(configfile):
    rawxmlfile = open(configfile).read()
    filename = xmlparser.parse(rawxmlfile,"FileName")
    rootfilenamelist = ", "
    for file in filename:
        if not file[:5]=="/eos/": rootfilenamelist += file+", "
    return rootfilenamelist[:-2]

def createcondortxt(jobname, jobnumber, jobdir):
    rootfiles = getoutputfilenames(jobname+"/xml/"+jobname+"_"+str(jobnumber)+".xml")
    additionalfiles = getinputfilenames(jobname+"/xml/"+jobname+"_"+str(jobnumber)+".xml")
    os.chdir(jobname+"/configs")
    condorfile = open(jobname+"_"+str(jobnumber)+".txt", 'w')
    outputdir = options.output
    if outputdir == "" : outputdir = jobdir+"/results"
    print >> condorfile, """universe = vanilla
Executable = %s/configs/%s_%d.sh
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
InitialDir = %s
Transfer_Input_Files = %s/configs/%s.tgz%s
Transfer_Output_Files = %s
Output = %s/logs/%s_%d.stdout
Error = %s/logs/%s_%d.stderr
Log = %s/logs/%s_%d.log
notify_user = ${LOGNAME}@FNAL.GOV
Arguments = %s 0
Queue 1""" %(jobname, jobname, jobnumber, outputdir, jobdir, jobname, additionalfiles, rootfiles.replace(".root","."+str(jobnumber)+".root"), jobdir, jobname, jobnumber, jobdir, jobname, jobnumber, jobdir, jobname, jobnumber, os.getcwd())
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
FILENAME=%s_%d.xml
cp %s/xml/${FILENAME} .
sed -i 's|FileName="/[^e][^o][^s].*/\(.*.root\)"|FileName="./\1"|' $FILENAME
mv ${WORKINGDIR}/*.root .
echo 'Running' >& $STATUSFILE
sframe_main %s_%d.xml
for filename in `/bin/ls *.root`; do
    newfilename=`echo $filename | sed 's|.root|.%d.root|'`
    mv $filename $newfilename
done
mv *.root $WORKINGDIR
echo 'Done' >& $STATUSFILE""" %(eosstatusdir, jobname, jobnumber, jobname, jobnumber, jobname, jobname, jobnumber, jobnumber)
    os.chmod(scriptname, 493) #493==755 in python chmod
    os.chdir("../..")

def begindatablock(mydatablock):
    return '<InputData Lumi="'+mydatablock.blocklumi+'" NEventsMax="'+mydatablock.neventsmax+'" Type="'+mydatablock.type+'" Version="'+mydatablock.version+'" Cacheable="'+mydatablock.cacheable+'">\n'

def enddatablack(mydatablock,indent):
    returnstring = indent+"  "+'<InputTree Name="AnalysisTree" />\n'
    if len(mydatablock.namelist)>1: returnstring += indent+"  "+'<OutputTree Name="AnalysisTree" />\n'
    returnstring += indent+"</InputData>\n"
    return returnstring

def makedatablocks(xmlfile,options):
    input=xmlfile
    datablocklist=[]
    while input.find("<InputData ") != -1:
        datablock = input[input.find("<InputData "):input.find("</InputData>")+12]
        version = xmlparser.parse(datablock,"Version")[0]
        vetoflag=0
        if options.veto.find(",") != -1: vetolist=options.veto.split(",")
        else: vetolist = [options.veto]
        if options.veto != "":
            for veto in vetolist:
                vetoflag = re.search(veto,version)
                if vetoflag: break
        filterflag=1
        if options.filter.find(",") != -1: filterlist=options.filter.split(",")
        else: filterlist = [options.filter]
        if options.filter != "":
            for filter in filterlist:
                filterflag = re.search(filter,version)
                if filterflag: break
        if filterflag and not vetoflag: datablocklist.append(datablock)
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
    #xmldump=open("BSFrame_config_dump.xml",'w');xmldump.write(input);xmldump.close()
    return input

def createxmlfiles(options):
    pdfmax = 1
    if options.pdf == "CT10": pdfmax = 52
    elif options.pdf == "cteq66": pdfmax = 44
    rawxmlfile = open(options.configxml).read()
    xmlfile = resolveentities(rawxmlfile)
    xmldatablocks = makedatablocks(xmlfile,options)
    totalfilelist = []
    datablocklist = []
    for xmldatablock in xmldatablocks:
        type = xmlparser.parse(xmldatablock,"Type")[0]
        version = xmlparser.parse(xmldatablock,"Version")[0]
        maxevents = xmlparser.parse(xmldatablock,"NEventsMax")[0]
        cacheable = xmlparser.parse(xmldatablock,"Cacheable")[0]
        filelist = xmlparser.parse(xmldatablock,"FileName")
        totalfilelist.extend(filelist)
        lumi = xmlparser.parse(xmldatablock,"Lumi")
        blocklumi = lumi.pop(0)
        namelist = xmlparser.parse(xmldatablock,"Name")
        mydatablock=datablock.datablock(blocklumi, filelist, lumi, namelist, type, version, maxevents, cacheable)
        datablocklist.append(mydatablock)
    if options.numjobs==0 and options.pdf=="": options.numjobs=len(totalfilelist)
    elif options.numjobs==0 and options.pdf!="": options.numjobs=1
    numjoblist=[]
    for i in range(options.numjobs): numjoblist.append(0)
    jobindex = 0
    for i in range(len(totalfilelist)):
        if jobindex>=options.numjobs: jobindex-=options.numjobs
        numjoblist[jobindex]+=1
        jobindex+=1
    for pdfindex in range(pdfmax):
        datablocknumber = 0
        blockindex = 0
        for jobnumber in range(options.numjobs):
            blockindex,datablocknumber = createxmlfile(xmlfile, jobnumber+1, datablocklist, datablocknumber, blockindex, numjoblist[jobnumber], options, pdfindex+1)

def createxmlfile(infile, jobnumber, datablocklist, datablocknumber, blockindex, numfiles, options, pdfindex):
    jobnumber = jobnumber+options.numjobs*(pdfindex-1)
    filename = options.jobname+"_"+str(jobnumber)+".xml"
    os.chdir(options.jobname+"/xml")
    if options.ttbargencut: infile = additem(infile, "ApplyMttbarGenCut", "True")
    if options.ttbargencut: infile = applyttbargencut(infile)
    if options.flavor != "": infile = applyflavorselection(infile, options.flavor)
    if options.jec != "": infile = applyjesystematic(infile, "JEC", options.jec)
    if options.jer != "": infile = applyjesystematic(infile, "JER", options.jer)
    if options.pileupfile != "": infile = changepileupfile(infile, options.pileupfile)
    if options.bjets != "": infile = applybjetsystematic(infile, options.bjets)
    if options.tjets != "": infile = applytjetsystematic(infile, options.tjets)
    if options.pdf != "": infile = applypdfsystematics(infile, options, pdfindex)
    frontend = infile[:infile.find("<InputData ")]
    indent = frontend[frontend.rfind("\n")+1:]
    backend = infile[infile.rfind("</InputData>")+13:]
    if len(datablocklist)>0:
        inputfilestring = begindatablock(datablocklist[datablocknumber])
        filelist=datablocklist[datablocknumber].filelist
        lumilist=datablocklist[datablocknumber].lumilist
    else:
        print "Error: No DataBlocks Found!\n"
        exit(3)

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
            inputfilestring += indent+"  "+'<In FileName="'+filelist[blockindex]+'" Lumi="'+lumilist[blockindex]+'"/>\n'
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
    return -1

def checklog(jobname, jobnumber):
    errorinfo = os.popen("egrep -i 'exit|break|exceed|error|traceback|aborted' "+jobname+"/logs/"+jobname+"_"+str(jobnumber)+".log").readline().strip("\n")
    errorline = errorinfo[errorinfo.find(")")+2:]
    return errorline

def checkstdout(jobname, jobnumber):
    errors = os.popen('egrep -i "exit|break|exceed|error|traceback|aborted|E R R O R|find tree AnalysisTree|fatal" '+jobname+"/logs/"+jobname+"_"+str(jobnumber)+".stdout").readlines()
    returnerror = ""
    if len(errors)>0:
        error = errors[0].strip("\n")
        if error.find(":") != -1: returnerror = error.split(":")[-1:][0]
        else: returnerror=error
    return returnerror

def getjobinfo(jobname,jobnumber,resubmitjobs):
    outputfiles = os.popen("grep Transfer_Output_Files "+jobname+"/configs/"+jobname+"_"+str(jobnumber)+".txt").readline().strip("\n")
    outputfiles = outputfiles.split(" ")[2:]
    outputdirectory = os.popen("grep InitialDir "+jobname+"/configs/"+jobname+"_"+str(jobnumber)+".txt | awk '{print $3}'").readline().strip("\n")
    jobinfo=""
    for file in outputfiles:
        file=file.strip(",")
        filepath = outputdirectory+"/"+file
        jobstatus = open(options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status").read().strip("\n")
        if os.path.isfile(filepath):
            if os.path.getsize(filepath)==0:
                if jobstatus=="Done":
                    jobinfo += " Root File "+file+" Empty!"
                    if resubmitjobs.count(jobnumber)<1: resubmitjobs.append(jobnumber)
            else:
                rootfile = ROOT.TFile.Open(filepath)
                try: iszombie=rootfile.IsZombie()
                except: iszombie=True
                if iszombie:
                    jobinfo += " Root File "+file+" is Zombie!"
                    if resubmitjobs.count(jobnumber)<1: resubmitjobs.append(jobnumber)
                elif rootfile.Get("AnalysisTree"):
                    analysistree = rootfile.Get("AnalysisTree")
                    jobinfo += " Root File "+file+" is Valid: "+str(int(analysistree.GetEntries()))+" Events."
                else:
                    hist = rootfile.Get("nprocessed")
                    jobinfo += " Root File "+file+" is Valid: "+str(int(hist.GetEntries()))+" Events."
                    if file.find("PostSelection") == -1: jobinfo += " Warning No AnalysisTree Found in "+file
                if not iszombie: rootfile.Close()
        else:
            jobinfo += " Output file "+file+" is not found!"
            if resubmitjobs.count(jobnumber)<1: resubmitjobs.append(jobnumber)
        if jobstatus=="Done":
            stdouterror = checkstdout(jobname, jobnumber)
            logerror = ""
            #logerror = checklog(jobname, jobnumber)
            if (stdouterror != "" or logerror != "") and resubmitjobs.count(jobnumber)<1: resubmitjobs.append(jobnumber)
            if stdouterror != "": jobinfo += " "+stdouterror
            if logerror != "": jobinfo += " "+logerror
    return jobinfo

if not options.create and options.submit=="" and options.kill=="" and not options.status:
    print "ERROR: Must either create, submit jobs, kill, or check the status of jobs"

workingdir=os.getcwd()
cmsswbase=os.getenv("CMSSW_BASE")
username=os.getenv("USER")
currentnode=os.getenv("HOST")
eosstatusdir="/eos/uscms/store/user/"+username+"/BSFrameStatus/"+options.jobname

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

    if not os.path.isdir(options.output) and options.output!="": os.makedirs(options.output)
    if not os.path.isdir("/eos/uscms/store/user/"+username+"/BSFrameStatus"): os.mkdir("/eos/uscms/store/user/"+username+"/BSFrameStatus")
    if not os.path.isdir(eosstatusdir): os.mkdir(eosstatusdir)
    else:
        shutil.rmtree(eosstatusdir)
        os.mkdir(eosstatusdir)

    print "Creating configuration files for task: "+options.jobname
    createxmlfiles(options)
    for xmlfile in os.popen("/bin/ls "+options.jobname+"/xml/*.xml").readlines():
        xmlfile = xmlfile.strip("\n")
        jobnumber = int(xmlfile[xmlfile.rfind("_")+1:xmlfile.rfind(".")])
        createcondortxt(options.jobname,jobnumber,workingdir+"/"+options.jobname)
        createcondorscript(options.jobname,jobnumber,eosstatusdir)
        os.system("echo 'Created' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

    if not options.notar:
        os.chdir(cmsswbase+"/..")
        tarball = options.jobname+".tgz"
        target = os.popen("echo ${CMSSW_BASE##*/}").readline().strip("\n")+"/"
        print "Creating tarball of "+target+" area."
        os.system("tar -czf "+tarball+" "+target+" --exclude-caches")
        os.system("mv "+tarball+" "+workingdir+"/"+options.jobname+"/configs")
        os.chdir(workingdir+"/"+options.jobname)
        os.system('echo "Signature: 8a477f597d28d172789f06886806bc55" >& CACHEDIR.TAG')
    os.chdir(workingdir)

if not os.path.isdir(options.jobname):
    print "ERROR: Job directory "+options.jobname+" does not exist!\nPlease create job with bsframe.py -c myconfig.xml --create."
    exit(5)

if options.retar:
    if options.create: print "There is no point in creating a task and then recreating the tarball."
    if options.notar: print "You are stupid!"
    os.chdir(workingdir+"/"+options.jobname)
    if os.path.isfile("CACHEDIR.TAG"): os.remove("CACHEDIR.TAG")
    os.chdir(cmsswbase+"/..")
    tarball = options.jobname+".tgz"
    target = os.popen("echo ${CMSSW_BASE##*/}").readline().strip("\n")+"/"
    print "Creating tarball of "+target+" area."
    os.system("tar -czf "+tarball+" "+target+" --exclude-caches")
    os.system("mv "+tarball+" "+workingdir+"/"+options.jobname+"/configs")
    os.chdir(workingdir+"/"+options.jobname)
    os.system('echo "Signature: 8a477f597d28d172789f06886806bc55" >& CACHEDIR.TAG')
    os.chdir(workingdir)

if options.kill!="":
    joblist=[]
    if options.kill=="all":
        options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip('\n'))
        joblist=range(1,options.numjobs+1)
    else: joblist=makejoblist(options.kill)
    condorstatus=os.popen("condor_q -submitter $USER").read()
    for jobnumber in joblist:
        print "Killing job number: %d" %(jobnumber)
        logfile = os.popen("/bin/ls -rt "+options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log | tail -1").readline().strip('\n')
        jobid = os.popen("grep submitted "+logfile+" | tail -1 | awk '{print $2}'").readline().strip("\n()").split(".")[0]
        if condorstatus.find(jobid) == -1:
            print "Error! Job "+jobid+" not found."
            continue
        subnode = getcondornode(jobid)
        if subnode == -1: print "Error! Job condor submission node not found."
        else:
            if subnode==currentnode: os.system("condor_rm "+jobid)
            else: os.system("rsh "+subnode+" condor_rm "+jobid)
            time.sleep(0.3)
            os.system("echo 'Killed' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

if options.clean!="":
    joblist=[]
    if options.clean=="all":
        options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip('\n'))
        joblist=range(1,options.numjobs+1)
    else: joblist=makejoblist(options.submit)
    print "Cleaning %d jobs" %(len(joblist))
    for jobnumber in joblist:
        if os.path.isfile(options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".log"): os.system("/bin/rm "+options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".log")
        if os.path.isfile(options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".stderr"): os.system("/bin/rm "+options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".stderr")
        if os.path.isfile(options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".stdout"): os.system("/bin/rm "+options.jobname+"/log/"+options.jobname+"_"+str(jobnumber)+".stdout")
        if os.path.isfile(options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status"): os.system("/bin/rm "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")
        if os.path.isfile(eosstatusdir+"/"+options.jobname+"_"+str(jobnumber)+".status"): os.system("/bin/rm "+eosstatusdir+"/"+options.jobname+"_"+str(jobnumber)+".status")
        resultsdir = os.popen("grep InitialDir "+options.jobname()+"/configs/"+options.jobname()+"_"+jobnumber+".txt | awk '{print $3}'").readline().strip('\n')
        rootfiles = getoutputfilenames(options.jobname+"/xml/"+options.jobname+"_"+str(jobnumber)+".xml")
        for rootfile in rootfiles:
            if os.path.isfile(resultsdir+"/"+rootfile): os.system("/bin/rm "+resultsdir+"/"+rootfile)

if options.submit!="":
    joblist=[]
    if options.submit=="all":
        options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip('\n'))
        joblist=range(1,options.numjobs+1)
    else: joblist=makejoblist(options.submit)
    print "Submitting %d jobs" %(len(joblist))
    for jobnumber in joblist:
        print "Submitting job number: %d" %(jobnumber)
        if not os.path.isfile(options.jobname+"/xml/"+options.jobname+"_"+str(jobnumber)+".xml"):
            print "Error: No configuration file for jobs number "+str(jobnumber)+"!"
            exit(6)
        subnum = int(os.popen("grep Arguments "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt | awk '{print $4}'").readline().strip('\n'))
        os.system("sed -i 's|Transfer_Output_Files = \(.*\)_"+str(subnum)+".root$|Transfer_Output_Files = \\1_"+str(subnum+1)+".root|' "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        os.system("sed -i 's|/configs "+str(subnum)+"$|/configs "+str(subnum+1)+"|' "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        if os.path.isfile(options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log"): os.system("/bin/rm "+options.jobname+"/logs/"+options.jobname+"_"+str(jobnumber)+".log")
        os.system("condor_submit "+options.jobname+"/configs/"+options.jobname+"_"+str(jobnumber)+".txt")
        os.system("echo 'Submitted' >& "+options.jobname+"/status/"+options.jobname+"_"+str(jobnumber)+".status")

jobstatuslist=[]
resubmitjobs=[]
if (options.status):
    print "Loading Root"
    import ROOT
    for statuslog in os.popen("/bin/ls "+eosstatusdir).readlines():
        statuslog = statuslog.strip("\n")
        if os.path.isfile(eosstatusdir+"/"+statuslog):
            eostimestamp = os.path.getmtime(eosstatusdir+"/"+statuslog)
            localtimestamp = os.path.getmtime(options.jobname+"/status/"+statuslog)
            if eostimestamp>localtimestamp: os.system("/bin/cp "+eosstatusdir+"/"+statuslog+" "+options.jobname+"/status/")
    print "\nJob Status Summary for Task: ",options.jobname
    print "================================================================================"
    print "Job Number         Status             Additional Information"
    print "--------------------------------------------------------------------------------"
    whitespace="                                                                                "
    options.numjobs=int(os.popen("/bin/ls "+options.jobname+"/xml/"+options.jobname+"_*.xml | wc -l").readline().strip("\n"))
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
