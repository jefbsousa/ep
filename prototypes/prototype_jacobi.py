#!/usr/bin/env python
# coding: utf-8

# https://www.quantstart.com/articles/Jacobi-Method-in-Python-and-NumPy

# https://www3.nd.edu/~zxu2/acms40390F12/Lec-7.3.pdf

# https://www.ucg.ac.me/skladiste/blog_10701/objava_23569/fajlovi/jacobi.pdf

# In[2]:


import numpy as np


# In[3]:


ITERATION_LIMIT = 1000
# initialize the matrix
A = np.array([[10., -1., 2., 0.],
 [-1., 11., -1., 3.],
 [2., -1., 10., -1.],
 [0.0, 3., -1., 8.]])
# initialize the RHS vector
b = np.array([6., 25., -11., 15.])


# In[4]:


A


# In[5]:


b


# In[7]:


print("System:" )

for i in range(A.shape[0]):
    row = ["{}*x{}" .format(A[i, j], j + 1) for j in range(A.shape[1])]
    print(" + ".join(row), "=", b[i])
print()


# In[25]:


np.dot(A[2, :2], x[:2])


# In[26]:


A[2, :2]


# In[27]:


x[:2]


# In[ ]:





# In[8]:


# np.zeros(len(b))
x = np.zeros_like (b)

for it_count in range(ITERATION_LIMIT):
    print("Current solution:" , x)
    x_new = np.zeros_like (x)
    
    for i in range(A.shape[0]):
        s1 = np.dot(A[i, :i], x[:i])
        s2 = np.dot(A[i, i + 1:], x[i + 1:])
        x_new[i] = (b[i] - s1 - s2) / A[i, i]
        
    if np.allclose(x, x_new, atol=1e-10, rtol=0.):
        break
        
    x = x_new
    


# In[9]:


# prints the system

print("Solution:" )
print(x)
error = np.dot(A, x) - b
print("Error:" )
print(error)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




