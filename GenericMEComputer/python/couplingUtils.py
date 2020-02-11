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


## ArXiv 1906.06949 eq 5.3, 5.4, 5.5
def _SMHLoops_F12(tau):
    return -2*tau-2*tau*(1-tau)*_SMHLoops_ftau(tau)
def _SMHLoops_F1(tau):
    return 2+3*tau*(1+(2-tau)*_SMHLoops_ftau(tau))
def _SMHLoops_ftau(tau):
    if tau >= 1:
        return (math.asin(math.sqrt(1/tau)))**2
    else:
        rad = math.sqrt(1-tau)
        return -0.25*(complex(math.log((1+rad)/(1-rad)), -math.pi))**2

def make_warsaw_smeftsim_MW(warsaw):
    ret = dict(p for p in warsaw.iteritems())
    zeropars = [    'cG',  'cGtil',  
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
                    'cuu1']
    for p in zeropars:
        if p not in ret: ret[p] = 0
    if 'LambdaSMEFT' not in ret: 
        ret['LambdaSMEFT'] = 1e3
    for p in ret.keys():
        if p not in zeropars and p != 'LambdaSMEFT':
            raise RuntimeError("Bogus parameter %s not in SMEFTsim Warsaw basis" % p)
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

def warsaw_smeftsim_to_smeftnlo(smeftsim):
    src = make_warsaw_smeftsim_MW(smeftsim)
    out = dict()
    out[('dim6','Lambda')] = src['LambdaSMEFT']
    out[('dim6','cpDC')  ] = src['cHDD']
    out[('dim6','cpWB')  ] = src['cHWB']
    out[('dim6','cdp')   ] = src['cHbox']
    out[('dim6','cp')    ] = src['cH']
    out[('dim6','cWWW')  ] = src['cW']
    out[('dim6','cG')    ] = src['cHWB']
    out[('dim6','cpG')   ] = src['cHG']
    out[('dim6','cpW')   ] = src['cHW']
    out[('dim6','cpBB')  ] = src['cHB']
    # I'm to lazy now to do the rest
    return out
    

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

def full_hc(hc):
    full = dict(p for p in hc.iteritems())
    zeropars = [ 'cosa', 'kHtt', 'kAtt', 'kHbb', 'kAbb', 'kHll', 'kAll', 'kHaa', 'kAaa', 'kHza', 'kAza', 'kHgg', 'kAgg', 'kHzz', 'kAzz', 'kHww', 'kAww', 'kHda', 'kHdz', 'kHdwR', 'kHdwI', 'kHcc', 'kAcc' ]
    for p in zeropars:
        if p not in full: full[p] = 0
    if 'Lambda' not in full:
        full['Lambda'] = 1e3
    if 'kSM' not in full:
        full['kSM'] = 1e3
    for k in full.iterkeys():
        if k not in zeropars and k not in ('Lambda', 'kSM'):
            raise RuntimeError("Bogus parameter %s not in HC basis" % k)
    return full

def hc2hig(hc, constants = constants_alphaScheme):
    ret = []
    all_constants = make_allconstants(constants)
    cw2, sw2, g2, gp2 = [ all_constants[x]**2 for x in ("cw","sw","g","gprime") ]

def hig2hc(hig, Lambda=1e3, alpha=0, constants = constants_alphaScheme):
    hig = make_hig(hig, constants = constants)

    cpval = max(abs(hig[k]) for k in 'tcAA tcZZ tcZA'.split())
    if cpval > 1e-6:
        if alpha == 0: 
            alpha = math.atan(cpval)
            print "ERROR, can't implement cp-violating operators with alpha = 0. Will set alpha = %g" % alpha

    all_constants = make_allconstants(constants)
    v, g, gp = [ all_constants[x]**2 for x in ("v","g","gprime") ]
    cw2, sw2, g2, gp2, e2 = [ all_constants[x]**2 for x in ("cw","sw","g","gprime","e") ]
    ca = math.cos(alpha); sa = math.sin(alpha)
    pi2 = math.pi**2

    ret = dict(Lambda=Lambda, cosa = ca,
            kSM = (1 + hig['dcZ'])/ca,
            kHaa = -72*pi2/(47*ca) * hig['cAA'],
            kAaa = - 3*pi2 / sa * hig['tcAA'] if sa else 0,
            kHgg = 12*pi2/ca * hig['cGG'],
            kAgg = - 8*pi2 / sa * hig['tcGG'] if sa else 0,
            kHzz = - (g2+gp2)*Lambda/(v*ca) * hig['cZZ'],
            kAzz = - (g2+gp2)*Lambda/(v*sa) * hig['tcZZ'] if sa else 0,
            kHww = - (  g2  )*Lambda/(v*ca) * (hig[ 'cZZ'] + 2*sw2*hig[ 'cZA'] + sw2*sw2*hig[ 'cAA']),
            kAww = - (  g2  )*Lambda/(v*sa) * (hig['tcZZ'] + 2*sw2*hig['tcZA'] + sw2*sw2*hig['tcAA'])  if sa else 0,
            kHdz = - (  g2  )*Lambda/(v*ca) * hig['cZBox'],
            kHda = - ( g*gp )*Lambda/((g2-gp2)*v*ca) * (
                                2*g2*hig['cZBox'] + (g2+gp2)*hig['cZZ'] - e2*hig['cAA'] - (g2-gp2)*hig['cZA']),
            kHdwR = - ( g2 )*Lambda/((g2-gp2)*v*ca) * (
                                g2*hig['cZBox'] + gp2*hig['cZZ'] - e2*sw2*hig['cAA'] - (g2-gp2)*sw2*hig['cZA']),
            kHdwI = 0., ## FIXME not sure what this is?
            kHza  = -144*pi2/(94*cw2-13)/ca * hig['cZA'],
            kAza  = -24*pi2/(8*cw2-5)/sa * hig['tcZA']  if sa else 0,
            )
    return full_hc(ret)


def make_hig(hig = {}, constants = constants_alphaScheme):
    """Set to zero missing coefficients, calculate dependent coefficients"""
    ## First, set to zero any missing coefficient
    ret = dict(p for p in hig.iteritems())
    for p in [ "cGG", "dcZ", "cAA", "cZA", "cZZ", "cZBox", 
               "tcGG", "tcAA", "tcZA", "tcZZ", 
               "dm" ]:
        if p not in ret: ret[p] = 0
    # now expand if necessary
    # YR4 arXiv 1610.07922 eq II.2.38 p297  (and checked against Rosetta)
    all_constants = make_allconstants(constants)
    cw, sw, mZ, g, gp, e = [ all_constants[x] for x in ("cw","sw","mZ","g","gprime","e") ]
    g2, gp2, e2 = g**2, gp**2, e**2
    exp = dict()
    exp['dcW'] = ret['dcZ'] + 4*ret['dm']
    exp['cWW'] = ret['cZZ']      + 2 * (sw**2) * ret['cZA']      + (sw**4) * ret['cAA'] 
    exp['tcWW'] = ret['tcZZ'] + 2 * (sw**2) * ret['tcZA'] + (sw**4) * ret['tcAA'] 
    exp['cWBox'] = (    g2 * ret['cZBox'] +    gp2   * ret['cZZ'] - e2*(sw**2) * ret['cAA'] - (g2-gp2) * (sw**2) * ret['cZA'] )/(g2-gp2)
    exp['cABox'] = (2 * g2 * ret['cZBox'] + (g2+gp2) * ret['cZZ'] - e2         * ret['cAA'] - (g2-gp2)           * ret['cZA'] )/(g2-gp2)
    for (k,v) in exp.iteritems():
        if k in ret:
            if abs(ret[k] - v) > 1e-9:
                raise RuntimeError("Parameter set already expanded: %s, old = %g new = %g diff = %g, all %s" % (k,ret[k],v,ret[k]-v,hig))
        else:
            ret[k] = v
    return ret

def make_free_higlike(hig = {}, constants = constants_alphaScheme):
    ## Set to zero any missing coefficient
    ret = dict(p for p in hig.iteritems())
    for p in [ "cGG", "dcZ", "cAA", "cZA", "cZZ", "cZBox", 
               "tcGG", "tcAA", "tcZA", "tcZZ", 
               "dm", "dcW" , "cWW", "tcWW", "cWBox", "cABox" ]:
        if p not in ret: ret[p] = 0
    return ret

def hig2jhu(hig, constants = constants_alphaScheme, SMLoop=False, Flip=False, noGamma=False):
    # JHU manual, eq 3
    all_constants = make_allconstants(constants)
    (mW, mZ, e, cw, sw, gs, mH) = [ all_constants[x] for x in ("mW","mZ","e","cw","sw","gs",'mH') ]
    mtop = 172.
    if SMLoop:
        SM_cAA_eff = - (_SMHLoops_F1((2*mW/mH)**2) + 3*(2/3.)**2*_SMHLoops_F12((2*mtop/mH)**2)) / (8 * math.pi**2)
        SM_cZA_eff = - 0.0592 ## FIXME analytic expression would be nice, but it's complicated
        if SMLoop == "flip":
            SM_cZA_eff = -SM_cZA_eff
    else:
        SM_cAA_eff = 0
        SM_cZA_eff = 0
    ###
    return dict(
            g1zz = 2.0  +     2        * hig['dcZ'],
            g2zz = -0.5 * (e/sw/cw)**2 * hig['cZZ'],
            l1zz =        (e/mZ/sw)**2 * hig['cZBox'], # this is k_{1}^{ZZ} / ( Lambda_{1}^{ZZ} )**2, in units of 1/GeV^2
            g4zz = -0.5 * (e/cw/sw)**2 * hig['tcZZ'],
            g1ww = 2.0 +             2 * hig['dcW'],
            g2ww = -0.5 * (e/sw)**2    * hig['cWW'],
            l1ww =        (e/mW/sw)**2 * hig['cWBox'],
            g4ww = -0.5 * (e/sw)**2    * hig['tcWW'],
            g2za = -0.5 * (e**2/sw/cw) *(hig['cZA'] + SM_cZA_eff) * (-1 if Flip else +1) * (0 if noGamma else 1),
            l1za = (e**2/sw/cw/mZ**2)  * hig['cABox']             * (-1 if Flip else +1) * (0 if noGamma else 1),
            g4za = -0.5 * (e**2/sw/cw) * hig['tcZA']              * (-1 if Flip else +1) * (0 if noGamma else 1),
            g2aa = -0.5 * (e**2)       *(hig['cAA'] + SM_cAA_eff) * (0 if noGamma else 1), 
            g4aa = -0.5 * (e**2)       * hig['tcAA']              * (0 if noGamma else 1),
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
        if abs(v) < 1e-8: continue
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


def byhand_hig_strings(scan, constants=constants_alphaScheme, SMLoop=False, Flip=False, noGamma=False):
    return [ "Name:%s %s" % (name, jhu2string(hig2jhu(point, constants=constants, SMLoop=SMLoop, Flip=Flip, noGamma=noGamma))) for (name,point) in scan ]


def _try_smeftsim(**kwargs):
    print sorted((k,v) for (k,v) in warsaw_smeftsim_MW_to_higgs(kwargs, constants_mWscheme).iteritems() if abs(v) > 1e-9)
def _try_smeftsim_tilde(**kwargs):
    scaleds = dict((k,v/0.24621965079413738**2) for (k,v) in kwargs.iteritems())
    print sorted((k,v) for (k,v) in warsaw_smeftsim_MW_to_higgs(scaleds, constants_mWscheme).iteritems() if abs(v) > 1e-9)
def _try_hel(**kwargs):
    print sorted((k,v) for (k,v) in hel2hig(kwargs).iteritems() if abs(v) > 1e-9)
