#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 15:40:59 2018

@author: christian
"""

import re


# function to clear up syntax of mathematica output
def repl_rate(pattern):
    pattern = pattern[0][0]+'('+pattern[0][1]+',1)'
    return(pattern)
def rpl_time_dep(pattern):
    pattern = '('+pattern[0][0]+')'
    return(pattern)
def clear_equation(eq):
    # cut off everything before equalty sign
    eq = re.sub('.+=','dydt(x) =',eq)
    # remove derivative expression
    eq = re.sub('-\(eta[0-9]\^\′\)\[t\]','',eq)
    # replace rate constants
    eq = re.sub('c[0-9]+',repl_rate,eq)
    eq = re.sub('(?<=/c\([0-9],)1','2',eq)
    # replace the time argument
    eq = re.sub('[0-9]+\[t\]',rpl_time_dep,eq)
    # get rit of rectangular bracket
    eq = re.sub('\[','(',eq)
    eq = re.sub('\]',')',eq)
    # replace capital by small l
    eq = re.sub('L','l',eq)
    # replace space by *
    eq = re.sub('(?<=[^=]) (?=[^=])','*',eq)
    return(eq)
def clear_grad(eq):
    # cut off everything before equalty sign
    eq = re.sub('.+=','grad(:,x) =',eq)
    # remove logarithm expression 
    eq = re.sub('[+-]Log\[.+c[0-d]\] m[0-9]\[t\]','',eq)
    eq = re.sub('[+-]Log\[.+c[0-d]\] \([^)]+\)','',eq)
    # replace the time argument
    eq = re.sub('[0-9]+\[t\]',rpl_time_dep,eq)
    # get rid of rectangular bracket
    eq = re.sub('\[','(',eq)
    eq = re.sub('\]',')',eq)
    # add colon as first argument
    eq = re.sub('m\((?=[0-9]+)','m(:,',eq)
    eq = re.sub('eta\((?=[0-9]+)','eta(:,',eq)
    # replace space by .*
    eq = re.sub('(?<=[^=]) (?=[^=])','.*',eq)
    # replace ^ by .^
    eq = re.sub('\^','.^',eq)
    return(eq)

file_folder = ''
file_name = 'variational_extended_gene_expression_central.txt'
file_path = file_folder+file_name
# load the file
file = open(file_path,'r')

# read text
content = file.read()
# produce list of lines
line_list = content.split('\n')


# create backward equations 
beq_ind_list = []
ident = 'beq'
for i in range(0,len(line_list)):
    if ident in line_list[i]:
        beq_ind_list.append(i)
# get the backward equations into a list
beq_list = []
for i in range(0,len(beq_ind_list)):
    tmp = line_list[beq_ind_list[i]+len(beq_ind_list)];
    tmp = clear_equation(tmp);
    tmp = re.sub('x',str(i+1),tmp)
    beq_list.append(tmp)
 
# create gradient equations
ceq_ind_list = []
ident = 'ceq'
for i in range(0,len(line_list)):
    if ident in line_list[i]:
        ceq_ind_list.append(i)
ceq_list = []
for i in range(0,len(ceq_ind_list)):
    tmp = line_list[ceq_ind_list[i]+len(ceq_ind_list)];
    tmp = clear_grad(tmp);
    tmp = re.sub('x',str(i+1),tmp)
    ceq_list.append(tmp)

# close input file file
file.close()

#create file for backward function
file_path = file_folder+'backward_equation.m'
file = open(file_path,'w+')
# print header stuff
file.write('function dydt = backward_equation(t,eta,alpha_t,alpha,m,c)'+'\n')
file.write('% The constraint equation in forward form'+'\n\n')
file.write('% perform interpolation'+'\n')
file.write('alpha = interp_new(alpha_t(2)-alpha_t(1),alpha,t-alpha_t(1));'+'\n')
file.write('m = interp_new(alpha_t(2)-alpha_t(1),m,t-alpha_t(1));'+'\n\n')
file.write('% evaluate the derivative'+'\n')
file.write('dydt = zeros('+str(len(beq_list))+',1);'+'\n')
# print equations to file
for beq in beq_list:
    file.write(beq+';\n')
# print end of function
file.write('\n')
file.write('end')
# close file
file.close()

# create file for forward function 
file_path = file_folder+'control_gradient.m'
file = open(file_path,'w+')
# print header stuff
file.write('function [grad] = control_gradient(eta,m)'+'\n')
file.write('% The constraint equation in forward form'+'\n\n')
file.write('% evalute the contribution to the gradient'+'\n')
file.write('grad = zeros(size(m,1),'+str(len(ceq_list))+');'+'\n');
# print equations to file
for ceq in ceq_list:
    file.write(ceq+';\n')
# print end of function
file.write('\n')
file.write('grad = -grad;'+'\n\n')
file.write('end')
# close file
file.close()

