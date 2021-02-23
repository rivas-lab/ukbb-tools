#!/bin/bash

ls ldsc/*log >log_list

bash ../07_LDSC/helpers/ldsc_rg_view.sh --list log_list

rm log_list
