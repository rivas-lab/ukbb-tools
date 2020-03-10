#!/bin/python
import json
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify

def check_credentials():
    # Debugging.  Sick of logging in and getting redirected
    #if True:
    #    return True
    access_token = session.get('access_token')
    if access_token is None:
        return False

    access_token = access_token[0]
    from urllib2 import Request, urlopen, URLError

    headers = {'Authorization': 'OAuth ' + access_token}
    req = Request('https://www.googleapis.com/oauth2/v1/userinfo',
                  None, headers)
    try:
        res = urlopen(req)
        data = json.loads(res.read())
        # todo log email address data['email']
        return True
    except URLError as e:
        if e.code == 401:
            # Unauthorized - bad token
            session.pop('access_token', None)
            return redirect(url_for('login_page'))
        return False

def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn


#@app.route('/mrp/<key>')
def mrp(key):
    #if not check_credentials():
    #    return redirect(url_for('login'))
    #db = get_db()

    with open(key + ".mcmc.posteriors", "r") as inFile:
        probabilities = []
        admixture_data = []
        variants = []
        for line in inFile.readlines():
            line = line.strip()
            var_info = line.split("\t")
            print('var_info')
            print(var_info)
            varid = var_info[0]
            print('varid')
            print(varid)
            var = var_info[4]
            print('var')
            print(var)
            variant = {}
            varidpage = varid.split(':')[0] + ':' + varid.split(':')[1]
            print('varidpage')
            print(varidpage)
            variant["variant"] = varidpage
            variant["varid"] = var
            for j in range(5, len(var_info)):
                variant["%s" % (j-4)] = var_info[j]
            variants.append(variant)
            print("variant",variant)
        print("variants",variants)
        admixture_data.append(variants)
    with open(key + ".mcmc.gene.posteriors", "r") as inFile:
        admixture_datagene = []
        genes = []
        for line in inFile.readlines():
            line = line.strip()
            gene_info = line.split("\t")
            geneid = gene_info[0]
            gene = {}
            gene["gene"] = geneid
            print(len(gene_info))
            for j in range(1, int((len(gene_info)-1)/3)+1):
                gene["%s" % (j-1)] = gene_info[j]
            genes.append(gene)
        admixture_datagene.append(genes)
    print("admixturedatagene",admixture_datagene)
    with open(key + ".lbf","r") as inFile:
        for line in inFile.readlines():
            line = line.strip()
            line = line.split()
            l10bf = "{0:.3g}".format(float(line[0])/numpy.log(10))
    with open(key + ".mcmc.bc", "r") as inFile:
        cluster_data = []
        headerarr = []
        inFiler = inFile.readlines()
        header = inFiler[0]
        header = header.split('\t')
        for j in range(1,len(header)):
            if j % 3 == 1:
#                headerarr.append(j)
                headeritem = header[j].decode('utf8').replace(" ","")
                headeritem = headeritem.rstrip(header[j][-3:]).upper()
                headerarr.append(headeritem)
        for line in inFiler[1:]:
            line = line.strip()
            cluster = []
            clusterl95 = []
            clusteru95 = []
            clust_info = line.split("\t")
            cluster_num = clust_info[0]
            for j in range(1, len(clust_info)):
                if j % 3 == 1:
                    tmp = []
                    headeritem = header[j].decode('utf8').replace(" ","")
                    headeritem = headeritem.rstrip(header[j][-3:]).upper()
                    tmp.append(headeritem.encode("ascii"))
                    print(headeritem)
                    tmp.append(float(clust_info[j]))
                if j % 3 == 2:
                    tmp.append(float(clust_info[j]))
                if j % 3 == 0:
                    tmp.append(float(clust_info[j]))
                    cluster.append(tmp)
            cluster_data.append([cluster_num, cluster])
    fdrf = open(key + ".fdr", "r").readlines()
    fdrnum = fdrf[0].rstrip().split()[0]
    fdr_data = []
    for line in fdrf[1:]:
        line = line.strip()
        line = line.split()
        fdr_data.append(line[0])
    #try:
    #    return render_template(
    #        'mrp.html',
    #        plot_data=admixture_data,
    #        gene_data=admixture_datagene,
    #        num_figs=cluster_data,
    #        log10_bf=l10bf,
    #        cluster_num=str(int(cluster_num)+1),
    #        fdr=fdrnum,
    #        fdr_data=fdr_data
    #        phenid_arr=headerarr
    #        )
    #except Exception as e:
    #    print('Failed: %s' % e)
    #    abort(404)

if __name__ == "__main__":
    mrp('ANGPTL7_test')
