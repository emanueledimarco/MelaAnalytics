import math

# follow sec 3.1 of arxiv 1709.06492
vckm_triv = [[1,0,0],[0,1,0],[0,0,1]]
constants_mWscheme = dict(mZ = 91.1876, mW = 80.365, Gf = 1.1663787e-5, mH = 125.09, alphas = 0.1181, vckm = vckm_triv)
constants_alphaScheme = dict(mZ = 91.1876, alpha = 1.0/127.950, Gf = 1.1663787e-5, mH = 125.09, alphas = 0.1181, vckm = vckm_triv)

def make_allconstants(constants):
    ret = dict(mZ = constants['mZ'], Gf = constants['Gf'], mH = constants['mH'], vckm = constants['vckm'])
    ret['v'] = 1.0/math.sqrt(math.sqrt(2) * ret['Gf'])
    if 'alpha' in constants:
        ret['alpha'] = constants['alpha']
        ret['e']  = math.sqrt( 4*math.pi*ret['alpha'] );
        ret['sw'] = math.sqrt(0.5*(1 - math.sqrt(1 - 4*math.pi*ret['alpha']/math.sqrt(2)/ret['Gf']/ret['mZ']**2)))
        ret['cw'] = math.sqrt(1 - ret['sw']**2)
        ret['mW'] = ret['mZ'] * ret['cw']
    elif 'mW' in constants:
        ret['mW'] = constants['mW']
        ret['cw'] = ret['mW']/ret['mZ'] 
        ret['sw'] = math.sqrt(1 - ret['cw']**2)
        ret['e']  = 2 * ret['mW'] * math.sqrt( math.sqrt(2) * ret['Gf'] ) * ret['sw'];
        ret['alpha'] = ret['e']**2/(4*math.pi)
    else: 
        raise RuntimeError("Constants must contain alpha or mW, while we have %s" % constants)
    ret['gprime'] = ret['e']/ret['cw']
    ret['g']      = ret['e']/ret['sw']
    if 'alphas' in constants:
        ret['alphas'] = constants['alphas']
        ret['gs'] = math.sqrt( 4*math.pi*constants['alphas'])
    elif 'gs' in constants:
        ret['gs'] = constants['gs']
        ret['alphas'] = constants['gs']**2/(4*math.pi)
    return ret

def make_warsaw_smeftsim_MW(warsaw):
    ret = dict(p for p in warsaw.iteritems())
    for p in [      'cG',  'cGtil',  
                    'cH',  
                    'cHB',  'cHBtil',  
                    'cHDD',  
                    'cHG',  'cHGtil',  
                    'cHW',  'cHWtil',  
                    'cHWB', 'cHWBtil',  
                    'cHbox',  
                    'cHd',  
                    'cHe',  
                    'cHl1',  
                    'cHl3',  
                    'cHq1',  
                    'cHq3',  
                    'cHu',  
                    'cHudAbs',  'cHudPh',  
                    'cW',  'cWtil',  
                    'cdBAbs',  'cdBPh',  
                    'cdGAbs',  'cdGPh',  
                    'cdHAbs',  'cdHPh',  
                    'cdWAbs',  'cdWPh',  
                    'cdd',  
                    'cdd1',  
                    'ceBAbs',  'ceBPh',  
                    'ceHAbs',  'ceHPh',  
                    'ceWAbs',  'ceWPh',  
                    'ced',  
                    'cee',  
                    'ceu',  
                    'cld',  
                    'cle',  
                    'cledqAbs',  'cledqPh',  
                    'clequ1Abs',  'clequ1Ph',  
                    'clequ3Abs',  'clequ3Ph',  
                    'cll',  
                    'cll1',  
                    'clq1',  
                    'clq3',  
                    'clu',  
                    'cqd1',  
                    'cqd8',  
                    'cqe',  
                    'cqq1',  
                    'cqq11',  
                    'cqq3',  
                    'cqq31',  
                    'cqu1',  
                    'cqu8',  
                    'cquqd1Abs',  'cquqd1Ph',  
                    'cquqd8Abs',  'cquqd8Ph',  
                    'cuBAbs',  'cuBPh',  
                    'cuGAbs',  'cuGPh',  
                    'cuHAbs',  'cuHPh',  
                    'cuWAbs',  'cuWPh',  
                    'cud1',  
                    'cud8',  
                    'cuu',  
                    'cuu1']:
        if p not in ret: ret[p] = 0
    if 'LambdaSMEFT' not in ret: 
        ret['LambdaSMEFT'] = 1e3
    return ret
def warsaw_smeftsim_MW_to_higgs(warsaw, constants=constants_mWscheme):
    src = make_warsaw_smeftsim_MW(warsaw)
    all_constants = make_allconstants(constants)
    # first do some stupid translation to Rosetta notation 
    scale = all_constants['v']**2 / src['LambdaSMEFT'] ** 2
    W = dict()
    for k, v in src.items():
        vs = v * scale
        if k == 'LambdaSMEFT': 
            continue
        elif k.endswith('Ph'):
            continue
        elif k.endswith('til'):
            k = 't' + k[:-3]
        elif k.endswith('Abs'):
            k = k[:-3]; phi = src[k+'Ph']
            vs = complex(v * math.cos(phi), v * math.sin(phi))
        W[k] = vs
    ## then translate as per Rosetta
    gw, gp, gs = [ all_constants[x] for x in ("g","gprime","gs") ]
    gw2, gp2, gs2 = gw**2, gp**2, gs**2
    ## some four-fermion elements (FIXME: not sure any of these are correct)
    cll1221 = W['cll'] # ?? 
    WBxH3l11 = W['cHl3'].real
    WBxH3l22 = W['cHl3'].real # should be identical in the U(3)^5 case ?
    M = dict()
    M['dm'] = ( -W['cHWB']*gw*gp - W['cHDD']*gw2/4.  + ( cll1221 - 2.*WBxH3l11 - 2.*WBxH3l22)*gp2/4.) / (gw2-gp2)
    M['dcZ'] = (W['cHbox'] - W['cHDD']/4. + 3./4.*cll1221 - 3./2.*WBxH3l11 - 3./2.*WBxH3l22) 
    M['cZBox'] = (-cll1221/2. + W['cHDD']/2.  + WBxH3l11 + WBxH3l22)/gw2
    M['cGG']  = (4./gs2)*W['cHG'] 
    M['cAA']  = 4.*(W['cHW']/gw2 + W['cHB']/gp2 - W['cHWB']/gw/gp)
    M['cZZ']  = 4.*( gw2*W['cHW'] + gp2*W['cHB'] + gw*gp*W['cHWB'] )/(gw2+gp2)**2
    M['cZA']  = (4.*W['cHW'] - 4.*W['cHB'] - 2.*(gw2-gp2)/(gw*gp)*W['cHWB'])/(gw2+gp2)
    M['tcGG']  = (4./gs2)*W['tcHG'] 
    M['tcAA']  = 4.*(W['tcHW']/gw2 + W['tcHB']/gp2 - W['tcHWB']/gw/gp)
    M['tcZZ']  = 4.*( gw2*W['tcHW'] + gp2*W['tcHB'] + gw*gp*W['tcHWB'] )/(gw2+gp2)**2
    M['tcZA']  = (4.*W['tcHW'] - 4.*W['tcHB'] - 2.*(gw2-gp2)/(gw*gp)*W['tcHWB'])/(gw2+gp2)
    # EWK boson couplings (not used by MELA)
    def f(T3,Q): # [eqn (A.4)]
        Acoeff = - gw*gp/(gw2-gp2)*W['cHWB']
        Zcoeff = (cll1221/4.- WBxH3l11/2.  - WBxH3l22/2. - W['cHDD']/4.)
        return Acoeff*Q + Zcoeff*(T3 + Q*gp2/(gw2-gp2))

    # W/Z chiral coupling deviations
    M['dgLwl'] = (WBxH3l11 + f(1./2.,0.) - f(-1./2.,-1.))
    M['dgLze'] = (-1./2.*WBxH3l11 - 1./2.*W['cHl1'].real + f(-1./2.,-1.))
    M['dgRze'] = (-1./2.*W['cHe'] + f(0.,-1.))
    M['dgLzu'] = (1./2.*W['cHq3'] - 1./2.*W['cHq1'] + f(1./2.,2./3.))
    M['dgLzd'] = -1./2.*(W['cHq3'] + W['cHq1']) + f(-1./2.,-1./3.) # for V_{CKM} = delta_{ij} ....
    M['dgRzu'] = (- 1./2.*W['cHu'] + f(0.,2./3.))
    M['dgRzd'] = (- 1./2.*W['cHd'] + f(0.,-1./3.))
    M['dgRwq'] = - 1./2.*W['cHud'].real
    # dependent but not handled by make_hig yet
    M['dgLzv'] = M['dgLwl'] + M['dgLze']
    M['dgLwq'] = M['dgLzu'] - M['dgLzd'] # for V_{CKM} = delta_{ij}

    # and complete the Higgs basis
    return make_hig(M, constants=constants)

def make_hel(hel):
    ret = dict(p for p in hel.iteritems())
    if "cW" in ret and "cWW" in ret and ret["cW"] != ret["cWW"]:
        raise RuntimeError("cW and cWW are the same operator")
    elif "cW" in ret:
        ret["cWW"] = ret["cW"]
    elif "cWW" in ret:
        ret["cW"] = ret["cWW"]
    else:
        ret["cW"] = 0; ret["cWW"] = 0
    for p in [  "cH", "cT", "c6", 
                "cu", "cd", "cl", 
                "cB", "cHW", "cHB", "cA", "cG", 
                "cHQ", "cpHQ", "cHu", "cHd", "cHud", "cHL", "cpHL", "cHe",
                "cuB", "cuW", "cuG", "cdB", "cdW", "cdG", "clB", "clW",
                "c3W", "c3G", "c2W", "c2B", "c2G",
                "tcHW", "tcHB", "tcG", "tcA", "tc3W", "tc3G" ]:
        if p not in ret: ret[p] = 0
    return ret

def hel2hig(hel, constants = constants_alphaScheme):
    # Use coefficients names as in the HEL UFO, LHC HXSWG 2019-04, CMS PAS HIG-19-005
    ## First, set to zero any missing coefficient
    src = make_hel(hel)
    ## Then translate 
    # YR4 arXiv 1610.07922 eq  II.2.9, II.2.21-23 p289+ and checked against Rosetta
    all_constants = make_allconstants(constants)
    sw2, g2, gp2 = [ all_constants[x]**2 for x in ("sw","g","gprime") ]
    cpHL22 = 0 # FIXME is this cpHL ?
    ret = dict()
    ret['dm'] = -gp2/(g2-gp2) * (src['cW'] + src['cB'] + src['c2W'] + src['c2B'] - g2/(2*gp2) * src['cT'] + 0.5*cpHL22)
    ret['dcZ'] = -0.5 * src['cH'] - 1.5*cpHL22
    ret['cGG'] = 16/(g2) * src['cG'] ## Note: it's g and not gs, since the gs is added in the lagrangian on top of it
    ret['cAA'] = 16/(g2) * src['cA']
    ret['cZZ'] = -4/(g2+gp2) * ( src['cHW'] +  gp2/g2 * src['cHB'] - 4*gp2/g2 * sw2 * src['cA'] )
    ret['cZBox'] = 2/g2 * ( src['cW'] + src['cHW'] + src['c2W'] + 
                             gp2/g2 * ( src['cB'] + src['cHB'] + src['c2B']) - 0.5*src['cT'] + 0.5*cpHL22 )
    ret['cZA'] = 2/g2 * ( src['cHB'] - src['cHW'] - 8 * sw2 * src['cA'] )
    ret['tcGG'] = 16/(g2) * src['tcG'] 
    ret['tcAA'] = 16/(g2) * src['tcA']
    ret['tcZZ'] = -4/(g2+gp2) * ( src['tcHW'] +  gp2/g2 * src['tcHB'] - 4*gp2/g2 * sw2 * src['tcA'] )
    ret['tcZA'] = 2/g2 * ( src['tcHB'] - src['tcHW'] - 8 * sw2 * src['tcA'] )
    # and complete the Higgs basis
    return make_hig(ret, constants=constants)

def make_hig(hig = {}, constants = constants_alphaScheme):
    """Set to zero missing coefficients, calculate dependent coefficients"""
    ## First, set to zero any missing coefficient
    ret = dict(p for p in hig.iteritems())
    for p in [ "cGG", "dcZ", "cAA", "cZA", "cZZ", "cZBox", 
               "tcGG", "tcAA", "tcZA", "tcZZ", 
               "dm" ]:
        if p not in ret: ret[p] = 0
    # now expand
    # YR4 arXiv 1610.07922 eq II.2.38 p297  (and checked against Rosetta)
    for k in "dcW" , "cWW", "tcWW", "cWBox", "cABox":
        if k in hig: raise RuntimeError("Parameter set already expanded: %s" % hig)
    all_constants = make_allconstants(constants)
    cw, sw, mZ, g, gp, e = [ all_constants[x] for x in ("cw","sw","mZ","g","gprime","e") ]
    g2, gp2, e2 = g**2, gp**2, e**2
    ret['dcW'] = ret['dcZ'] + 4*ret['dm']
    ret['cWW']      = ret['cZZ']      + 2 * (sw**2) * ret['cZA']      + (sw**4) * ret['cAA'] 
    ret['tcWW'] = ret['tcZZ'] + 2 * (sw**2) * ret['tcZA'] + (sw**4) * ret['tcAA'] 
    ret['cWBox'] = (    g2 * ret['cZBox'] +    gp2   * ret['cZZ'] - e2*(sw**2) * ret['cAA'] - (g2-gp2) * (sw**2) * ret['cZA'] )/(g2-gp2)
    ret['cABox'] = (2 * g2 * ret['cZBox'] + (g2+gp2) * ret['cZZ'] - e2         * ret['cAA'] - (g2-gp2)           * ret['cZA'] )/(g2-gp2)
    return ret

def make_free_higlike(hig = {}, constants = constants_alphaScheme):
    ## Set to zero any missing coefficient
    ret = dict(p for p in hig.iteritems())
    for p in [ "cGG", "dcZ", "cAA", "cZA", "cZZ", "cZBox", 
               "tcGG", "tcAA", "tcZA", "tcZZ", 
               "dm", "dcW" , "cWW", "tcWW", "cWBox", "cABox" ]:
        if p not in ret: ret[p] = 0
    return ret

def hig2jhu(hig, constants = constants_alphaScheme):
    # JHU manual, eq 3
    all_constants = make_allconstants(constants)
    (mW, mZ, e, cw, sw, gs) = [ all_constants[x] for x in ("mW","mZ","e","cw","sw","gs") ]
    return dict(
            g1zz = 2.0  +     2        * hig['dcZ'],
            g2zz = -0.5 * (e/sw/cw)**2 * hig['cZZ'],
            l1zz =        (e/mZ/sw)**2 * hig['cZBox'], # this is k_{1}^{ZZ} / ( Lambda_{1}^{ZZ} )**2, in units of 1/GeV^2
            g4zz = -0.5 * (e/cw/sw)**2 * hig['tcZZ'],
            g1ww = 2.0 +             2 * hig['dcW'],
            g2ww = -0.5 * (e/sw)**2    * hig['cWW'],
            l1ww =        (e/mW/sw)**2 * hig['cWBox'],
            g4ww = -0.5 * (e/sw)**2    * hig['tcWW'],
            g2za = -0.5 * (e**2/sw/cw) * hig['cZA'],
            l1za = (e**2/sw/cw/mZ**2)  * hig['cABox'],
            g4za = -0.5 * (e**2/sw/cw) * hig['tcZA'],
            g2aa = -0.5 * (e**2)       * hig['cAA'],
            g4aa = -0.5 * (e**2)       * hig['tcAA'],
            g2gg = -0.5 * (gs**2)      * hig['cGG'],
            g4gg = -0.5 * (gs**2)      * hig['tcGG'],
            )

def jhu2mela(jhu):
    # note: Lambda1's are all fixed to 10 TeV
    return dict(ghz1 = jhu['g1zz'], ghw1 = jhu['g1ww'],
                ghz2 = jhu['g2zz'], ghw2 = jhu['g2ww'], ghzgs2 = jhu['g2za'], ghgsgs2 = jhu['g2aa'],
                ghz4 = jhu['g4zz'], ghw4 = jhu['g4ww'], ghzgs4 = jhu['g4za'], ghgsgs4 = jhu['g4aa'],
                ghz1_prime2 = jhu['l1zz']*1e8, 
                ghw1_prime2 = jhu['l1ww']*1e8, 
                ghzgs1_prime2 = jhu['l1za']*1e8)

def mela2string(mela):
    ret = ["separateWWZZcouplings=1"] 
    for k,v in mela.iteritems():
        if v == 0: continue
        if k[0] == "g": 
            if isinstance(v, complex):
                ret.append("%s=%g,%g" % (k,v.real,v.imag))
            else:
                ret.append("%s=%g,0" % (k,v))
        else          : 
            ret.append("%s=%g"   % (k,v))
    return "Couplings:" + ";".join(ret)

def jhu2string(jhu): 
    return mela2string(jhu2mela(jhu))

def make_scan(couplings,quadratic=False):
    ret = [ ("SM",dict()) ]
    for (c,v) in couplings.iteritems():
        ret.append((c+"_up", dict([(c,+v)])))
        ret.append((c+"_dn", dict([(c,-v)])))
    if quadratic:
        for i,(c1,v1) in enumerate(couplings):
            for (c2,v2) in enumerate(couplings[(i+1):]):
                ret.append(("%s_%s_up" % (c1,c2), dict([(c1,v1),(c2,v2)])))
    return ret

def hig_scan_strings(couplings, constants=constants_alphaScheme, quadratic=False):
    scan = make_scan(couplings, quadratic=quadratic)
    return [ "Name:%s %s" % (name, jhu2string(hig2jhu(make_hig(point), constants=constants))) for (name,point) in scan ]

def hel_scan_strings(couplings, constants=constants_alphaScheme, quadratic=False):
    scan = make_scan(couplings, quadratic=quadratic)
    return [ "Name:%s %s" % (name, jhu2string(hig2jhu(hel2hig(point), constants=constants))) for (name,point) in scan ]

def warsaw_scan_strings(couplings, constants=constants_mWscheme, quadratic=False):
    scan = make_scan(couplings, quadratic=quadratic)
    return [ "Name:%s %s" % (name, jhu2string(hig2jhu(warsaw_smeftsim_MW_to_higgs(point, constants=constants), constants=constants))) for (name,point) in scan ]


def byhand_hig_strings(scan, constants=constants_alphaScheme):
    return [ "Name:%s %s" % (name, jhu2string(hig2jhu(point, constants=constants))) for (name,point) in scan ]



