#/usr/bin/python

from shutil import *
import os

import math
import shutil
import os
import subprocess
import getpass
import time
#import commands
import sys
import tempfile

queueName = "sht.q"

def GetEnvironmentVariable():
    env_var_list = os.environ
    
    # Now filter out problematic variables (the one from shellshock)
    if 'BASH_FUNC_module()' in env_var_list :
        #env_var_list.remove('BASH_FUNC_module()')
        #del os.environ['BASH_FUNC_module()']
        #os.unsetenv('BASH_FUNC_module()')
        env_var_list.pop('BASH_FUNC_module()')
    #env_var_list.remove('_')
    #env_var_list.remove('PS1')
    return env_var_list
    
def GetJobsInQueue():
    uname= getpass.getuser()
    stat = "qstat -u " + uname
    f=os.popen(stat)
    
    n=0
    for i in f.readlines():
      #  if i.find("mpi-hd.mpg.de")>-1:
        n=n+1
    int2return = n-2
    if int2return<0:
        return 0
    
    return int2return


def GetJobsInQueueName(nameFrag):

    uname= getpass.getuser()
    stat = "qstat -u " + uname 
    f=os.popen(stat)
    
    n=0
    for i in f.readlines():
        if i.find(nameFrag)>-1:
            n=n+1
   
    print(n,"Jobs containing string",nameFrag,"in queue")

    return n
 
def SubmitJob(cmd,jobname, extraopt):

    username = os.environ.get('USER')
    tmpdir2use = '/usr/tmp/{0}_py/'.format(username)

    if (not os.path.exists(tmpdir2use)):
        try:
            os.makedirs(tmpdir2use)
        except OSError:
            print("Can't create the directory [{0}] ==> Already exist or permission problem.".format(tmpdir2use))
    
    cwd = os.getcwd()
    #print(cwd)
    
    #filename="./hess-XXXX.sh"
    #print(filename)
    #libpath = os.environ.get("LD_LIBRARY_PATH")
    #path = os.environ.get("PATH")
    #pypath = os.environ.get("PYTHONPATH")

    fileprefix = username+'-'
    filesuffix = '.sh'
    f = tempfile.NamedTemporaryFile(dir=tmpdir2use,prefix=fileprefix,suffix=filesuffix,delete=False,mode='w')
    filename = f.name
    f.write("#!/bin/bash")
    f.write('\n')
    f.write('echo -------------------------------------\n')
    f.write('date\n')
    f.write('echo "USER    : $USER" \n')
    f.write('echo "JOB_ID  : $JOB_ID" \n')
    f.write('echo "JOB_NAME: $JOB_NAME" \n')
    f.write('echo "HOSTNAME: $HOSTNAME" \n')
    f.write('echo -------------------------------------\n')
    f.write("ulimit -c 1\n") # prevent big core dumps!
    #f.write("export LD_LIBRARY_PATH="+libpath+" \n")
    #f.write("export PATH="+path+" \n")
    #f.write("export PYTHONPATH="+pypath+" \n")
    #f.write('### Export environment variable (Brute Force Approach)\n')
    env_var_list = GetEnvironmentVariable()
    #print(env_var_list)
    for myenv in env_var_list:
        #print('export {0}=\"{1}\" \n'.format(myenv,os.environ.get(myenv)) )
        f.write(('export {0}=\"{1}\" \n'.format(myenv,os.environ.get(myenv))) )

    f.write(('cd {0} \n'.format(cwd)))
    #print(cmd+'\n')
    f.write(('{0}'.format(cmd)))
    f.write(('\n'))
    
    f.close()

    defaultopt = " -o "+cwd+"/ -e "+cwd+"/" ##+" -j y" 
    optfromuser = '';
    if (len(extraopt)==0) :
        optfromuser = defaultopt
    else:
        optfromuser = extraopt
    callCmd = ["qsub","-notify","-N" ,jobname, "-o", "$QUEUE_STDOUT", "-e", "$QUEUE_STDERR",filename]
    
    #qcmd = "qsub -q "+ queueName +" -notify " + " -N " + jobname + " " + optfromuser + " "+ filename
    qcmd = "qsub -P short "+" -notify " + " -N " + jobname + " " + optfromuser + " "+ filename #for short queue
    #qcmd = "qsub "+" -notify " + " -N " + jobname + " " + optfromuser + " "+ filename 
    #cmd = "qsub "+ filename

    print(qcmd)
    #subprocess.call(callCmd)
    #os.popen(cmd)

    #qsub_status = commands.getstatusoutput(qcmd)
    qsub_status = subprocess.call(qcmd,shell=True)
    print(qsub_status)

    #time.sleep(1)

  #  os.remove(filename) # Not needed it's in the /usr/tmp directory now

def SubmitJobWhenReady(cmd,jobname,maxjobs,extraopt=""):
    
    nj = GetJobsInQueue()
    print(nj,"out of",maxjobs,"allowed jobs currently submitted")

    while True:
        nj = GetJobsInQueue()

        if (nj<maxjobs):
            print("Submitting Job",jobname)
            SubmitJob(cmd,jobname,extraopt)
            break
        else:
            time.sleep(60)

#SubmitJobWhenReady("Test123","job",1)

if __name__ == '__main__':
    jobname=os.environ.get('USER')+'_AwesomeJob'
    maxjobs = 600
    input_arg_list = sys.argv[1:]
    #    print(input_arg_list)
    #    print(len(input_arg_list) )
    if (len(input_arg_list)==0):
        print("Usage : {0} MyCommandToSubmit (Be carefull with special caracters, better put everything in one string)".format(sys.argv[0]))
    else:
        mycmdtosubmit = str()
        for arginput in input_arg_list:
            mycmdtosubmit+=arginput+' '
        print("I'm going to submit this command : [ {0}]".format(mycmdtosubmit))
        SubmitJobWhenReady(mycmdtosubmit,jobname,maxjobs)
