#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Deleting files that we do not need, files generated from a previous run
import os
import glob

extension_list=["*inp*","*out","*ener","*rest*","*Hess*","*REST*","*.xyz","*.pdb"]
for name in extension_list:

    filelist=glob.glob(name)
    for file in filelist:
        os.remove(file)

