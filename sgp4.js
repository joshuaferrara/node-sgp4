var SGP4 = {
    deg2rad: Math.PI / 180.0,
    rad2deg: 57.2957795,
    twopi: 2.0 * Math.PI,
    fmod: function (x, y) {
      return x - Math.floor(x/y) * y;
    },
    dpper: function(satrec, inclo, init, ep, inclp, nodep, argpp, mp, opsmode) {
        'use strict';
        
        var e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, t, xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos, zns, zes, znl, zel;
        e3 = satrec.e3;
        ee2 = satrec.ee2;
        peo = satrec.peo;
        pgho = satrec.pgho;
        pho = satrec.pho;
        pinco = satrec.pinco;
        plo = satrec.plo;
        se2 = satrec.se2;
        se3 = satrec.se3;
        sgh2 = satrec.sgh2;
        sgh3 = satrec.sgh3;
        sgh4 = satrec.sgh4;
        sh2 = satrec.sh2;
        sh3 = satrec.sh3;
        si2 = satrec.si2;
        si3 = satrec.si3;
        sl2 = satrec.sl2;
        sl3 = satrec.sl3;
        sl4 = satrec.sl4;
        t = satrec.t;
        xgh2 = satrec.xgh2;
        xgh3 = satrec.xgh3;
        xgh4 = satrec.xgh4;
        xh2 = satrec.xh2;
        xh3 = satrec.xh3;
        xi2 = satrec.xi2;
        xi3 = satrec.xi3;
        xl2 = satrec.xl2;
        xl3 = satrec.xl3;
        xl4 = satrec.xl4;
        zmol = satrec.zmol;
        zmos = satrec.zmos;
        
        zns = 1.19459e-5;
        zes = 0.01675;
        znl = 1.5835218e-4;
        zel = 0.05490;
        
        var zm = zmos + zns * t;
        
        if (init === 'y') {
            zm = zmos;
        }
            
        var zf = zm + 2.0 * zes * Math.sin(zm);
        var sinzf = Math.sin(zf);
        var f2 = 0.5 * sinzf * sinzf - 0.25;
        var f3    = -0.5 * sinzf * Math.cos(zf);
        var ses   = se2* f2 + se3 * f3;
        var sis   = si2 * f2 + si3 * f3;
        var sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
        var sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
        var shs   = sh2 * f2 + sh3 * f3;
        zm    = zmol + znl * t;
        if (init === 'y') {
            zm = zmol;
        }
        zf    = zm + 2.0 * zel * Math.sin(zm);
        sinzf = Math.sin(zf);
        f2    =  0.5 * sinzf * sinzf - 0.25;
        f3    = -0.5 * sinzf * Math.cos(zf);
        var sel   = ee2 * f2 + e3 * f3;
        var sil   = xi2 * f2 + xi3 * f3;
        var sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
        var sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
        var shll  = xh2 * f2 + xh3 * f3;
        var pe    = ses + sel;
        var pinc  = sis + sil;
        var pl    = sls + sll;
        var pgh   = sghs + sghl;
        var ph    = shs + shll;
        
        if (init === 'n') {
            pe = pe - peo;
            pinc  = pinc - pinco;
            pl    = pl - plo;
            pgh   = pgh - pgho;
            ph    = ph - pho;
            inclp = inclp + pinc;
            ep    = ep + pe;
            var sinip = Math.sin(inclp);
            var cosip = Math.cos(inclp);
            
            if (inclp >= 0.2) {
                ph /= sinip;
                pgh -= cosip * ph;
                argpp += pgh;
                nodep += ph;
                mp += pl;
            } else {
                var sinop = Math.sin(nodep);
                var cosop = Math.cos(nodep);
                var alfdp = sinip * sinop;
                var betdp = sinip * cosop;
                var dalf = ph * cosop + pinc * cosip * sinop;
                var dbet = -ph * sinop + pinc * cosip * cosop;
                alfdp = alfdp + dalf;
                betdp = betdp + dbet;
                nodep = SGP4.fmod(nodep, SGP4.twopi);
                if (nodep < 0.0 && opsmode === 'a') {
                    nodep = nodep + SGP4.twopi;
                }
                var xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep;
                var xnoh = nodep;
                nodep = Math.atan2(alfdp, betdp);
                if (nodep < 0.0 && opsmode === 'a') {
                    nodep = nodep + SGP4.twopi;
                }
                if (Math.abs(xnoh - nodep) > Math.PI) {
                    if (nodep < xnoh) {
                        nodep = nodep + SGP4.twopi;
                    } else {
                        nodep = nodep - SGP4.twopi;
                    }
                }
                mp += pl;
                argpp = xls - mp - cosip * nodep;
            }
        }
        
        return {
            ep: ep,
            inclp: inclp,
            nodep: nodep,
            argpp: argpp,
            mp: mp
        };
    },
    dscom: function(epoch, ep, argpp, tc, inclp, nodep,  np, e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, xgh2, xgh3, xgh4,  xh2, xh3, xi2, xi3, xl2, xl3, xl4, zmol, zmos) {
        'use strict';
        
        var a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, betasq, cc,  ctem, stem, x1, x2, x3, x4, x5, x6, x7, x8, xnodce, xnoi, zcosg, zsing, zcosgl, zsingl, zcosh, zsinh, zcoshl, zsinhl, zcosi, zsini, zcosil, zsinil, zx, zy, ss1,  ss2,  ss3,  ss4,  ss5,  ss6,  ss7, sz1,  sz2,  sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33, s1, s2, s3, s4, s5, s6, s7, z1, z2, z3, z11,z12,z13, z21,z22,z23, z31,z32,z33;
        
        var zes     =  0.01675;
        var zel     =  0.05490;
        var c1ss    =  2.9864797e-6;
        var c1l     =  4.7968065e-7;
        var zsinis  =  0.39785416;
        var zcosis  =  0.91744867;
        var zcosgs  =  0.1945905;
        var zsings  = -0.98088458;
         
        var nm     = np;
        var em     = ep;
        var snodm  = Math.sin(nodep);
        var cnodm  = Math.cos(nodep);
        var sinomm = Math.sin(argpp);
        var cosomm = Math.cos(argpp);
        var sinim  = Math.sin(inclp);
        var cosim  = Math.cos(inclp);
        var emsq   = em * em;
        betasq = 1.0 - emsq;
        var rtemsq = Math.sqrt(betasq);
        
        peo    = 0.0;
        pinco  = 0.0;
        plo    = 0.0;
        pgho   = 0.0;
        pho    = 0.0;
        var day    = epoch + 18261.5 + tc / 1440.0;
        xnodce = (4.5236020 - 9.2422029e-4 * day) % SGP4.twopi;
        stem   = Math.sin(xnodce);
        ctem   = Math.cos(xnodce);
        zcosil = 0.91375164 - 0.03568096 * ctem;
        zsinil = Math.sqrt(1.0 - zcosil * zcosil);
        zsinhl = 0.089683511 * stem / zsinil;
        zcoshl = Math.sqrt(1.0 - zsinhl * zsinhl);
        var gam    = 5.8351514 + 0.0019443680 * day;
        zx     = 0.39785416 * stem / zsinil;
        zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
        zx     = Math.atan2(zx, zy);
        zx     = gam + zx - xnodce;
        zcosgl = Math.cos(zx);
        zsingl = Math.sin(zx);
        
        zcosg = zcosgs;
        zsing = zsings;
        zcosi = zcosis;
        zsini = zsinis;
        zcosh = cnodm;
        zsinh = snodm;
        cc    = c1ss;
        xnoi  = 1.0 / nm;
        
        for (var lsflg = 1; lsflg <= 2; lsflg++) {
            a1  =   zcosg * zcosh + zsing * zcosi * zsinh;
            a3  =  -zsing * zcosh + zcosg * zcosi * zsinh;
            a7  =  -zcosg * zsinh + zsing * zcosi * zcosh;
            a8  =   zsing * zsini;
            a9  =   zsing * zsinh + zcosg * zcosi * zcosh;
            a10 =   zcosg * zsini;
            a2  =   cosim * a7 + sinim * a8;
            a4  =   cosim * a9 + sinim * a10;
            a5  =  -sinim * a7 + cosim * a8;
            a6  =  -sinim * a9 + cosim * a10;
    
            x1  =  a1 * cosomm + a2 * sinomm;
            x2  =  a3 * cosomm + a4 * sinomm;
            x3  = -a1 * sinomm + a2 * cosomm;
            x4  = -a3 * sinomm + a4 * cosomm;
            x5  =  a5 * sinomm;
            x6  =  a6 * sinomm;
            x7  =  a5 * cosomm;
            x8  =  a6 * cosomm;
    
            z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
            z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
            z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
            z1  =  3.0 *  (a1 * a1 + a2 * a2) + z31 * emsq;
            z2  =  6.0 *  (a1 * a3 + a2 * a4) + z32 * emsq;
            z3  =  3.0 *  (a3 * a3 + a4 * a4) + z33 * emsq;
            z11 = -6.0 * a1 * a5 + emsq *  (-24.0 * x1 * x7-6.0 * x3 * x5);
            z12 = -6.0 * (a1 * a6 + a3 * a5) + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
            z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
            z21 =  6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
            z22 =  6.0 *  (a4 * a5 + a2 * a6) + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
            z23 =  6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
            z1  = z1 + z1 + betasq * z31;
            z2  = z2 + z2 + betasq * z32;
            z3  = z3 + z3 + betasq * z33;
            s3  = cc * xnoi;
            s2  = -0.5 * s3 / rtemsq;
            s4  = s3 * rtemsq;
            s1  = -15.0 * em * s4;
            s5  = x1 * x3 + x2 * x4;
            s6  = x2 * x3 + x1 * x4;
            s7  = x2 * x4 - x1 * x3;
            
            if (lsflg === 1) {
                ss1   = s1;
                ss2   = s2;
                ss3   = s3;
                ss4   = s4;
                ss5   = s5;
                ss6   = s6;
                ss7   = s7;
                sz1   = z1;
                sz2   = z2;
                sz3   = z3;
                sz11  = z11;
                sz12  = z12;
                sz13  = z13;
                sz21  = z21;
                sz22  = z22;
                sz23  = z23;
                sz31  = z31;
                sz32  = z32;
                sz33  = z33;
                zcosg = zcosgl;
                zsing = zsingl;
                zcosi = zcosil;
                zsini = zsinil;
                zcosh = zcoshl * cnodm + zsinhl * snodm;
                zsinh = snodm * zcoshl - cnodm * zsinhl;
                cc    = c1l;
            }
        }
        zmol = (4.7199672 + 0.22997150  * day - gam) % SGP4.twopi;
        zmos = (6.2565837 + 0.017201977 * day) % SGP4.twopi;
        
        se2  =   2.0 * ss1 * ss6;
        se3  =   2.0 * ss1 * ss7;
        si2  =   2.0 * ss2 * sz12;
        si3  =   2.0 * ss2 * (sz13 - sz11);
        sl2  =  -2.0 * ss3 * sz2;
        sl3  =  -2.0 * ss3 * (sz3 - sz1);
        sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * zes;
        sgh2 =   2.0 * ss4 * sz32;
        sgh3 =   2.0 * ss4 * (sz33 - sz31);
        sgh4 = -18.0 * ss4 * zes;
        sh2  =  -2.0 * ss2 * sz22;
        sh3  =  -2.0 * ss2 * (sz23 - sz21);
        
        ee2  =   2.0 * s1 * s6;
        e3   =   2.0 * s1 * s7;
        xi2  =   2.0 * s2 * z12;
        xi3  =   2.0 * s2 * (z13 - z11);
        xl2  =  -2.0 * s3 * z2;
        xl3  =  -2.0 * s3 * (z3 - z1);
        xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel;
        xgh2 =   2.0 * s4 * z32;
        xgh3 =   2.0 * s4 * (z33 - z31);
        xgh4 = -18.0 * s4 * zel;
        xh2  =  -2.0 * s2 * z22;
        xh3  =  -2.0 * s2 * (z23 - z21);
        
        return {
            snodm: snodm,
            cnodm: cnodm,
            sinim: sinim,
            cosim: cosim,
            sinomm: sinomm,
            cosomm: cosomm,
            day: day,
            e3: e3,
            ee2: ee2,
            em: em,
            emsq: emsq,
            gam: gam,
            peo: peo,
            pgho: pgho,
            pho: pho,
            pinco: pinco,
            plo: plo,
            rtemsq: rtemsq,
            se2: se2,
            se3: se3,
            sgh2: sgh2,
            sgh3: sgh3,
            sgh4: sgh4,
            sh2: sh2,
            sh3: sh3,
            si2: si2,
            si3: si3,
            sl2: sl2,
            sl3: sl3,
            sl4: sl4,
            s1: s1,
            s2: s2,
            s3: s3,
            s4: s4,
            s5: s5,
            s6: s6,
            s7: s7,
            ss1: ss1,
            ss2: ss2,
            ss3: ss3,
            ss4: ss4,
            ss5: ss5,
            ss6: ss6,
            ss7: ss7,
            sz1: sz1,
            sz2: sz2,
            sz3: sz3,
            sz11: sz11,
            sz12: sz12,
            sz13: sz13,
            sz21: sz21,
            sz22: sz22,
            sz23: sz23,
            sz31: sz31,
            sz32: sz32,
            sz33: sz33,
            xgh2: xgh2,
            xgh3: xgh3,
            xgh4: xgh4,
            xh2: xh2,
            xh3: xh3,
            xi2: xi2,
            xi3: xi3,
            xl2: xl2,
            xl3: xl3,
            xl4: xl4,
            nm: nm,
            z1: z1,
            z2: z2,
            z3: z3,
            z11: z11,
            z12: z12,
            z13: z13,
            z21: z21,
            z22: z22,
            z23: z23,
            z31: z31,
            z32: z32,
            z33: z33,
            zmol: zmol,
            zmos: zmos
        };
    },
    dsinit: function(whichconst, cosim, emsq, argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, t, tc, gsto, mo, mdot, no, nodeo, nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33, ecco, eccsq, em, argpm, inclm, mm, nm, nodem, irez, atime, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, didt, dmdt, dnodt, domdt, del1, del2, del3, xfact, xlamo, xli, xni) {
        'use strict';
        
        var f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543, g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533, sini2, temp, temp1, theta, xno2, ainv2,  aonv,   cosisq, eoc;
        
        var q22    = 1.7891679e-6;
        var q31    = 2.1460748e-6;
        var q33    = 2.2123015e-7;
        var root22 = 1.7891679e-6;
        var root44 = 7.3636953e-9;
        var root54 = 2.1765803e-9;
        var rptim  = 4.37526908801129966e-3; // equates to 7.29211514668855e-5 rad/sec
        var root32 = 3.7393792e-7;
        var root52 = 1.1428639e-7;
        var x2o3   = 2.0 / 3.0;
        var znl    = 1.5835218e-4;
        var zns    = 1.19459e-5;
        
        var xke = whichconst.xke;
        
        irez = 0;
        if (0.0034906585 < nm < 0.0052359877) {
            irez = 1;
        }
        if (8.26e-3 <= nm <= 9.24e-3 && em >= 0.5) {
            irez = 2;
        }
            
        var ses  =  ss1 * zns * ss5;
        var sis  =  ss2 * zns * (sz11 + sz13);
        var sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
        var sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
        var shs  = -zns * ss2 * (sz21 + sz23);
        // sgp4fix for 180 deg incl
        if (inclm < 5.2359877e-2 || inclm > Math.PI - 5.2359877e-2) {
            shs = 0.0;
        }
        if (sinim !== 0.0) {
            shs = shs / sinim;
        }
        var sgs  = sghs - cosim * shs;
        
        dedt = ses + s1 * znl * s5;
        didt = sis + s2 * znl * (z11 + z13);
        dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
        var sghl = s4 * znl * (z31 + z33 - 6.0);
        var shll = -znl * s2 * (z21 + z23);
        // sgp4fix for 180 deg incl
        if (inclm < 5.2359877e-2 || inclm > Math.PI - 5.2359877e-2) {
            shll = 0.0;
        }
        domdt = sgs + sghl;
        dnodt = shs;
        
        if (sinim !== 0.0) {
            domdt = domdt - cosim / sinim * shll;
            dnodt = dnodt + shll / sinim;
        }
        
        var dndt   = 0.0;
        theta  = (gsto + tc * rptim) % SGP4.twopi;
        em     = em + dedt * t;
        inclm  = inclm + didt * t;
        argpm  = argpm + domdt * t;
        nodem  = nodem + dnodt * t;
        mm     = mm + dmdt * t;
        
        if (irez !== 0) {
            aonv = Math.pow(nm / xke, x2o3);
            if (irez === 2) {
                cosisq = cosim * cosim;
                var emo    = em;
                em     = ecco;
                var emsqo  = emsq;
                emsq   = eccsq;
                eoc    = em * emsq;
                g201   = -0.306 - (em - 0.64) * 0.440;
    
                if (em <= 0.65) {
                    g211 =    3.616  -  13.2470 * em +  16.2900 * emsq;
                    g310 =  -19.302  + 117.3900 * em - 228.4190 * emsq +  156.5910 * eoc;
                    g322 =  -18.9068 + 109.7927 * em - 214.6334 * emsq +  146.5816 * eoc;
                    g410 =  -41.122  + 242.6940 * em - 471.0940 * emsq +  313.9530 * eoc;
                    g422 = -146.407  + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
                    g520 = -532.114  + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
                } else {
                    g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
                    g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
                    g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
                    g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
                    g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
                    if (em > 0.715) {
                        g520 =-5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc;
                    } else {
                        g520 = 1464.74 -  4664.75 * em +  3763.64 * emsq;
                    }
                }
    
                if (em < 0.7) {
                    g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21  * eoc;
                    g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                    g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4  * eoc;
                } else {
                    g533 =-37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
                    g521 =-51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                    g532 =-40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
                }
    
                sini2=  sinim * sinim;
                f220 =  0.75 * (1.0 + 2.0 * cosim+cosisq);
                f221 =  1.5 * sini2;
                f321 =  1.875 * sinim  *  (1.0 - 2.0 * cosim - 3.0 * cosisq);
                f322 = -1.875 * sinim  *  (1.0 + 2.0 * cosim - 3.0 * cosisq);
                f441 = 35.0 * sini2 * f220;
                f442 = 39.3750 * sini2 * sini2;
                f522 =  9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim- 5.0 * cosisq) + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq) );
                f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq) + 6.56250012 * (1.0+2.0 * cosim - 3.0 * cosisq));
                f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim+cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
                f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim+cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
                xno2  =  nm * nm;
                ainv2 =  aonv * aonv;
                temp1 =  3.0 * xno2 * ainv2;
                temp  =  temp1 * root22;
                d2201 =  temp * f220 * g201;
                d2211 =  temp * f221 * g211;
                temp1 =  temp1 * aonv;
                temp  =  temp1 * root32;
                d3210 =  temp * f321 * g310;
                d3222 =  temp * f322 * g322;
                temp1 =  temp1 * aonv;
                temp  =  2.0 * temp1 * root44;
                d4410 =  temp * f441 * g410;
                d4422 =  temp * f442 * g422;
                temp1 =  temp1 * aonv;
                temp  =  temp1 * root52;
                d5220 =  temp * f522 * g520;
                d5232 =  temp * f523 * g532;
                temp  =  2.0 * temp1 * root54;
                d5421 =  temp * f542 * g521;
                d5433 =  temp * f543 * g533;
                xlamo =  (mo + nodeo + nodeo-theta - theta) % SGP4.twopi;
                xfact =  mdot + dmdt + 2.0 * (nodedot + dnodt - rptim) - no;
                em    = emo;
                emsq  = emsqo;
            }
            if (irez === 1) {
                g200  = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
                g310  = 1.0 + 2.0 * emsq;
                g300  = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
                f220  = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                f311  = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
                f330  = 1.0 + cosim;
                f330  = 1.875 * f330 * f330 * f330;
                del1  = 3.0 * nm * nm * aonv * aonv;
                del2  = 2.0 * del1 * f220 * g200 * q22;
                del3  = 3.0 * del1 * f330 * g300 * q33 * aonv;
                del1  = del1 * f311 * g310 * q31 * aonv;
                xlamo = (mo + nodeo + argpo - theta) % SGP4.twopi;
                xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no;
            }
            xli   = xlamo;
            xni   = no;
            atime = 0.0;
            nm    = no + dndt;
        }
        return {
            em: em,
            argpm: argpm,
            inclm: inclm,
            mm: mm,
            nm: nm,
            nodem: nodem,
            irez: irez,
            atime: atime,
            d2201: d2201,
            d2211: d2211,
            d3210: d3210,
            d3222: d3222,
            d4410: d4410,
            d4422: d4422,
            d5220: d5220,
            d5232: d5232,
            d5421: d5421,
            d5433: d5433,
            dedt: dedt,
            didt: didt,
            dmdt: dmdt,
            dndt: dndt,
            dnodt: dnodt,
            domdt: domdt,
            del1: del1,
            del2: del2,
            del3: del3,
            xfact: xfact,
            xlamo: xlamo,
            xli: xli,
            xni: xni
        };
    },
    dspace: function(irez, d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, dedt, del1, del2, del3, didt, dmdt, dnodt, domdt, argpo, argpdot, t, tc, gsto, xfact, xlamo, no, atime, em, argpm, inclm, xli, mm, xni, nodem, nm) {
        'use strict';
        
        var delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi;
        
        var fasx2 = 0.13130908;
        var fasx4 = 2.8843198;
        var fasx6 = 0.37448087;
        var g22   = 5.7686396;
        var g32   = 0.95240898;
        var g44   = 1.8014998;
        var g52   = 1.0508330;
        var g54   = 4.4108898;
        var rptim = 4.37526908801129966e-3; // equates to 7.29211514668855e-5 rad/sec
        var stepp =    720.0;
        var stepn =   -720.0;
        var step2 = 259200.0;
        
        var dndt   = 0.0;
        theta  = (gsto + tc * rptim) % SGP4.twopi;
        em     = em + dedt * t;
    
        inclm  = inclm + didt * t;
        argpm  = argpm + domdt * t;
        nodem  = nodem + dnodt * t;
        mm     = mm + dmdt * t;
        
        ft    = 0.0;
        
        if (irez !== 0) {
            // sgp4fix streamline check
            if (atime === 0.0 || t * atime <= 0.0 || Math.abs(t) < Math.abs(atime)) {
                atime  = 0.0;
                xni    = no;
                xli    = xlamo;
            }
    
            // sgp4fix move check outside loop
            if (t > 0.0) {
                delt = stepp;
            } else {
                delt = stepn;
            }
    
            var iretn = 381; // added for do loop
            var iret  =   0; // added for loop
            while (iretn === 381) {
                //  ------------------- dot terms calculated -------------
                //  ----------- near - synchronous resonance terms -------
                if (irez !== 2) {
                    xndt  = del1 * Math.sin(xli - fasx2) + del2 * Math.sin(2.0 * (xli - fasx4)) + del3 * Math.sin(3.0 * (xli - fasx6));
                    xldot = xni + xfact;
                    xnddt = del1 * Math.cos(xli - fasx2) + 2.0 * del2 * Math.cos(2.0 * (xli - fasx4)) + 3.0 * del3 * Math.cos(3.0 * (xli - fasx6));
                    xnddt = xnddt * xldot;
                } else {
                    // --------- near - half-day resonance terms --------
                    xomi  = argpo + argpdot * atime;
                    x2omi = xomi + xomi;
                    x2li  = xli + xli;
                    xndt  = (d2201 * Math.sin(x2omi + xli - g22) + d2211 * Math.sin(xli - g22) + d3210 * Math.sin(xomi + xli - g32)  + d3222 * Math.sin(-xomi + xli - g32)+ d4410 * Math.sin(x2omi + x2li - g44)+ d4422 * Math.sin(x2li - g44) + d5220 * Math.sin(xomi + xli - g52)  + d5232 * Math.sin(-xomi + xli - g52)+ d5421 * Math.sin(xomi + x2li - g54) + d5433 * Math.sin(-xomi + x2li - g54));
                    xldot = xni + xfact;
                    xnddt = (d2201 * Math.cos(x2omi + xli - g22) + d2211 * Math.cos(xli - g22) + d3210 * Math.cos(xomi + xli - g32) + d3222 * Math.cos(-xomi + xli - g32) + d5220 * Math.cos(xomi + xli - g52) + d5232 * Math.cos(-xomi + xli - g52) + 2.0 * (d4410 * Math.cos(x2omi + x2li - g44) + d4422 * Math.cos(x2li - g44) + d5421 * Math.cos(xomi + x2li - g54) + d5433 * Math.cos(-xomi + x2li - g54)));
                    xnddt = xnddt * xldot;
                }
                //  ----------------------- integrator -------------------
                //  sgp4fix move end checks to end of routine
                if (Math.abs(t - atime) >= stepp) {
                    iret  = 0;
                    iretn = 381;
                } else {
                    ft    = t - atime;
                    iretn = 0;
                }
    
                if (iretn === 381) {
                    xli   = xli + xldot * delt + xndt * step2;
                    xni   = xni + xndt * delt + xnddt * step2;
                    atime = atime + delt;
                }
            }
            nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
            xl = xli + xldot * ft + xndt * ft * ft * 0.5;
            if (irez !== 1) {
                mm   = xl - 2.0 * nodem + 2.0 * theta;
                dndt = nm - no;
            } else {
                mm   = xl - nodem - argpm + theta;
                dndt = nm - no;
            }
             
            nm = no + dndt;
        }
        return {
            atime: atime,
            em: em,
            argpm: argpm,
            inclm: inclm,
            xli: xli,
            mm: mm,
            xni: xni,
            nodem: nodem,
            dndt: dndt,
            nm: nm
        };
    },
    initl: function(satn, whichconst, ecco, epoch, inclo, no, method, opsmode) {
        'use strict';
        
        var ak, d1,  del,  adel, po, gsto;
        
        var tumin = whichconst.tumin;
        var mu = whichconst.mu;
        var radiusearthkm = whichconst.radiusearthkm;
        var xke = whichconst.xke;
        var j2 = whichconst.j2;
        var j3 = whichconst.j3;
        var j4 = whichconst.j4;
        var j3oj2 = whichconst.j3oj2;
        
        var x2o3   = 2.0 / 3.0;
        
        var eccsq  = ecco * ecco;
        var omeosq = 1.0 - eccsq;
        var rteosq = Math.sqrt(omeosq);
        var cosio  = Math.cos(inclo);
        var cosio2 = cosio * cosio;
        
        ak    = Math.pow(xke / no, x2o3);
        d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
        var del_  = d1 / (ak * ak);
        adel  = ak * (1.0 - del_ * del_ - del_ *
                 (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
        del_  = d1/(adel * adel);
        no    = no / (1.0 + del_);
    
        var ao    = Math.pow(xke / no, x2o3);
        var sinio = Math.sin(inclo);
        po    = ao * omeosq;
        var con42 = 1.0 - 5.0 * cosio2;
        var con41 = -con42-cosio2-cosio2;
        var ainv  = 1.0 / ao;
        var posq  = po * po;
        var rp    = ao * (1.0 - ecco);
        method = 'n';
        
        if (opsmode === 'a') {
            //  sgp4fix use old way of finding gst
            //  count integer number of days from 0 jan 1970
            var ts70  = epoch - 7305.0;
            var ds70 = Math.floor(ts70 + 1.0e-8); // 1.0;
            var tfrac = ts70 - ds70;
            //  find greenwich location at epoch
            var c1    = 1.72027916940703639e-2;
            var thgr70= 1.7321343856509374;
            var fk5r  = 5.07551419432269442e-15;
            var c1p2p = c1 + SGP4.twopi;
            gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % SGP4.twopi;
            if (gsto < 0.0) {
                gsto = gsto + SGP4.twopi;
            }
        } else {
            gsto = SGP4.gstime(epoch + 2433281.5);
        }
        return {
            no: no,
            method: method,
            ainv: ainv,
            ao: ao,
            con41: con41,
            con42: con42,
            cosio: cosio,
            cosio2: cosio2,
            eccsq: eccsq,
            omeosq: omeosq,
            posq: posq,
            rp: rp,
            rteosq: rteosq,
            sinio: sinio ,
            gsto: gsto
        };
    },
    sgp4init: function(whichconst, opsmode, satn, epoch, xbstar, xecco, xargpo, xinclo, xmo, xno, xnodeo, satrec) {
        'use strict';
        
        var cnodm, snodm, cosim, sinim, cosomm, sinomm, cc1sq, cc2, cc3, coef, coef1, cosio4, day, dndt, em, emsq, eeta, etasq, gam, argpm, nodem, inclm, mm, nm, perige, pinvsq, psisq, qzms24, rtemsq, s1, s2, s3, s4, s5, s6, s7, sfour, ss1,ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33, tc, temp, temp1, temp2, temp3, temp4, tsi, xpidot, xhdot1, z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33;
        
        temp4 = 1.5e-12;
        satrec.isimp   = 0;   satrec.method = 'n'; satrec.aycof    = 0.0;
        satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0;
        satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0;
        satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0;
        satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0;
        satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0;
        satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0;
        satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot  = 0.0;
        satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf   = 0.0;
        
        satrec.irez  = 0;   satrec.d2201 = 0.0; satrec.d2211 = 0.0;
        satrec.d3210 = 0.0; satrec.d3222 = 0.0; satrec.d4410 = 0.0;
        satrec.d4422 = 0.0; satrec.d5220 = 0.0; satrec.d5232 = 0.0;
        satrec.d5421 = 0.0; satrec.d5433 = 0.0; satrec.dedt  = 0.0;
        satrec.del1  = 0.0; satrec.del2  = 0.0; satrec.del3  = 0.0;
        satrec.didt  = 0.0; satrec.dmdt  = 0.0; satrec.dnodt = 0.0;
        satrec.domdt = 0.0; satrec.e3    = 0.0; satrec.ee2   = 0.0;
        satrec.peo   = 0.0; satrec.pgho  = 0.0; satrec.pho   = 0.0;
        satrec.pinco = 0.0; satrec.plo   = 0.0; satrec.se2   = 0.0;
        satrec.se3   = 0.0; satrec.sgh2  = 0.0; satrec.sgh3  = 0.0;
        satrec.sgh4  = 0.0; satrec.sh2   = 0.0; satrec.sh3   = 0.0;
        satrec.si2   = 0.0; satrec.si3   = 0.0; satrec.sl2   = 0.0;
        satrec.sl3   = 0.0; satrec.sl4   = 0.0; satrec.gsto  = 0.0;
        satrec.xfact = 0.0; satrec.xgh2  = 0.0; satrec.xgh3  = 0.0;
        satrec.xgh4  = 0.0; satrec.xh2   = 0.0; satrec.xh3   = 0.0;
        satrec.xi2   = 0.0; satrec.xi3   = 0.0; satrec.xl2   = 0.0;
        satrec.xl3   = 0.0; satrec.xl4   = 0.0; satrec.xlamo = 0.0;
        satrec.zmol  = 0.0; satrec.zmos  = 0.0; satrec.atime = 0.0;
        satrec.xli   = 0.0; satrec.xni   = 0.0;
        
        satrec.bstar   = xbstar;
        satrec.ecco    = xecco;
        satrec.argpo   = xargpo;
        satrec.inclo   = xinclo;
        satrec.mo	    = xmo;
        satrec.no	    = xno;
        satrec.nodeo   = xnodeo;
        
        satrec.operationmode = opsmode;
        
        var tumin = whichconst.tumin;
        var mu = whichconst.mu;
        var radiusearthkm = whichconst.radiusearthkm;
        var xke = whichconst.xke;
        var j2 = whichconst.j2;
        var j3 = whichconst.j3;
        var j4 = whichconst.j4;
        var j3oj2 = whichconst.j3oj2;
        
        var ss = 78.0 / radiusearthkm + 1.0;
        var qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
        var qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
        var x2o3 =  2.0 / 3.0;
        
        satrec.init = 'y';
        satrec.t = 0.0;
        
        var method, ainv, ao, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio;
        satrec.no, method, ainv, ao, satrec.con41, con42, cosio, cosio2, eccsq, omeosq, posq, rp, rteosq, sinio, satrec.gsto;
        var initResults = SGP4.initl(satn, whichconst, satrec.ecco, epoch, satrec.inclo, satrec.no, satrec.method, satrec.operationmode);

        satrec.no       = initResults.no;
        method      = initResults.method;
        ainv        = initResults.ainv;
        ao          = initResults.ao;
        satrec.con41    = initResults.con41;
        con42       = initResults.con42;
        cosio       = initResults.cosio;
        cosio2      = initResults.cosio2;
        eccsq       = initResults.eccsq;
        omeosq      = initResults.omeosq;
        posq        = initResults.posq;
        rp          = initResults.rp;
        rteosq      = initResults.rteosq;
        sinio       = initResults.sinio;
        satrec.gsto     = initResults.gsto;
        
        satrec.error = 0;
        
        if (omeosq >= 0.0 || satrec.no >= 0.0) {
            satrec.isimp = 0;
            if (rp < 220.0 / radiusearthkm + 1.0) {
                satrec.isimp = 1;
            }
            sfour  = ss;
            qzms24 = qzms2t;
            perige = (rp - 1.0) * radiusearthkm;

            //for perigees below 156 km, s and qoms2t are altered -
            if (perige < 156.0) {
                sfour = perige - 78.0;
                if (perige < 98.0) {
                    sfour = 20.0;
                }
                //  sgp4fix use multiply for speed instead of pow
                var qzms24temp =  (120.0 - sfour) / radiusearthkm;
                qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                sfour  = sfour / radiusearthkm + 1.0;
            }
    
            pinvsq = 1.0 / posq;
    
            tsi  = 1.0 / (ao - sfour);
            satrec.eta  = ao * satrec.ecco * tsi;
            etasq = satrec.eta * satrec.eta;
            eeta  = satrec.ecco * satrec.eta;
            psisq = Math.abs(1.0 - etasq);
            coef  = qzms24 * Math.pow(tsi, 4.0);
            coef1 = coef / Math.pow(psisq, 3.5);
            cc2   = coef1 * satrec.no * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.375 * j2 * tsi / psisq * satrec.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
            satrec.cc1   = satrec.bstar * cc2;
            cc3   = 0.0;
            if (satrec.ecco > 1.0e-4) {
                cc3 = -2.0 * coef * tsi * j3oj2 * satrec.no * sinio / satrec.ecco;
            }
            satrec.x1mth2 = 1.0 - cosio2;
            satrec.cc4 = 2.0 * satrec.no * coef1 * ao * omeosq * (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecco * (0.5 + 2.0 * etasq) - j2 * tsi / (ao * psisq) * (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * Math.cos(2.0 * satrec.argpo)));
            satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
            cosio4 = cosio2 * cosio2;
            temp1  = 1.5 * j2 * pinvsq * satrec.no;
            temp2  = 0.5 * temp1 * j2 * pinvsq;
            temp3  = -0.46875 * j4 * pinvsq * pinvsq * satrec.no;
            satrec.mdot = satrec.no + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
            satrec.argpdot = (-0.5 * temp1 * con42 + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4) + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4));
            xhdot1 = -temp1 * cosio;
            satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
            xpidot =  satrec.argpdot+ satrec.nodedot;
            satrec.omgcof = satrec.bstar * cc3 * Math.cos(satrec.argpo);
            satrec.xmcof = 0.0;
            if (satrec.ecco > 1.0e-4) {
                satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta;
            }
            satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1;
            satrec.t2cof   = 1.5 * satrec.cc1;
            //  sgp4fix for divide by zero with xinco = 180 deg
            if (Math.abs(cosio+1.0) > 1.5e-12) {
                satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
            } else {
                satrec.xlcof = -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
            }
            satrec.aycof   = -0.5 * j3oj2 * sinio;
            // sgp4fix use multiply for speed instead of pow
            var delmotemp = 1.0 + satrec.eta * Math.cos(satrec.mo);
            satrec.delmo   = delmotemp * delmotemp * delmotemp;
            satrec.sinmao  = Math.sin(satrec.mo);
            satrec.x7thm1  = 7.0 * cosio2 - 1.0;

            // --------------- deep space initialization -------------
            if (2*Math.PI / satrec.no >= 225.0) {
                satrec.method = 'd';
                satrec.isimp  = 1;
                tc    =  0.0;
                inclm = satrec.inclo;
        
                //var z33, z32, z31, z23, z22, z21, z13, z12, z11, z3, z2, z1, nm, sz33, sz32, sz31, sz23, sz22, sz21, sz13, sz12, sz11, sz3, sz2, sz1, ss7, ss6, ss5, ss4, ss3, ss2, ss1, s7, s6, s5, s4, s3, s2, s1, rtemsq, gam, emsq, em, day, cosomm, sinomm, cosim, sinim, cnodm, snodm;
                
                var dscomResults = SGP4.dscom(epoch, satrec.ecco, satrec.argpo, tc, satrec.inclo, satrec.nodeo, satrec.no, satrec.e3, satrec.ee2, satrec.peo, satrec.pgho, satrec.pho, satrec.pinco, satrec.plo, satrec.se2, satrec.se3, satrec.sgh2, satrec.sgh3, satrec.sgh4, satrec.sh2, satrec.sh3, satrec.si2, satrec.si3, satrec.sl2, satrec.sl3, satrec.sl4, satrec.xgh2, satrec.xgh3, satrec.xgh4, satrec.xh2, satrec.xh3, satrec.xi2, satrec.xi3, satrec.xl2, satrec.xl3, satrec.xl4, satrec.zmol, satrec.zmos);
                
                snodm = dscomResults.snodm;
                cnodm = dscomResults.cnodm;
                sinim = dscomResults.sinim;
                cosim = dscomResults.cosim;
                sinomm = dscomResults.sinomm;
                cosomm = dscomResults.cosomm;
                day = dscomResults.day;
                satrec.e3 = dscomResults.e3;
                satrec.ee2 = dscomResults.ee2;
                em = dscomResults.em;
                emsq = dscomResults.emsq;
                gam = dscomResults.gam;
                satrec.peo = dscomResults.peo;
                satrec.pgho = dscomResults.pgho;
                satrec.pho = dscomResults.pho;
                satrec.pinco = dscomResults.pinco;
                satrec.plo = dscomResults.plo;
                rtemsq = dscomResults.rtemsq;
                satrec.se2 = dscomResults.se2;
                satrec.se3 = dscomResults.se3;
                satrec.sgh2 = dscomResults.sgh2;
                satrec.sgh3 = dscomResults.sgh3;
                satrec.sgh4 = dscomResults.sgh4;
                satrec.sh2 = dscomResults.sh2;
                satrec.sh3 = dscomResults.sh3;
                satrec.si2 = dscomResults.si2;
                satrec.si3 = dscomResults.si3;
                satrec.sl2 = dscomResults.sl2;
                satrec.sl3 = dscomResults.sl3;
                satrec.sl4 = dscomResults.sl4;
                s1 = dscomResults.s1;
                s2 = dscomResults.s2;
                s3 = dscomResults.s3;
                s4 = dscomResults.s4;
                s5 = dscomResults.s5;
                s6 = dscomResults.s6;
                s7 = dscomResults.s7;
                ss1 = dscomResults.ss1;
                ss2 = dscomResults.ss2;
                ss3 = dscomResults.ss3;
                ss4 = dscomResults.ss4;
                ss5 = dscomResults.ss5;
                ss6 = dscomResults.ss6;
                ss7 = dscomResults.ss7;
                sz1 = dscomResults.sz1;
                sz2 = dscomResults.sz2;
                sz3 = dscomResults.sz3;
                sz11 = dscomResults.sz11;
                sz12 = dscomResults.sz12;
                sz13 = dscomResults.sz13;
                sz21 = dscomResults.sz21;
                sz22 = dscomResults.sz22;
                sz23 = dscomResults.sz23;
                sz31 = dscomResults.sz31;
                sz32 = dscomResults.sz32;
                sz33 = dscomResults.sz33;
                satrec.xgh2 = dscomResults.xgh2;
                satrec.xgh3 = dscomResults.xgh3;
                satrec.xgh4 = dscomResults.xgh4;
                satrec.xh2 = dscomResults.xh2;
                satrec.xh3 = dscomResults.xh3;
                satrec.xi2 = dscomResults.xi2;
                satrec.xi3 = dscomResults.xi3;
                satrec.xl2 = dscomResults.xl2;
                satrec.xl3 = dscomResults.xl3;
                satrec.xl4 = dscomResults.xl4;
                nm = dscomResults.nm;
                z1 = dscomResults.z1;
                z2 = dscomResults.z2;
                z3 = dscomResults.z3;
                z11 = dscomResults.z11;
                z12 = dscomResults.z12;
                z13 = dscomResults.z13;
                z21 = dscomResults.z21;
                z22 = dscomResults.z22;
                z23 = dscomResults.z23;
                z31 = dscomResults.z31;
                z32 = dscomResults.z32;
                z33 = dscomResults.z33;
                satrec.zmol = dscomResults.zmol;
                satrec.zmos = dscomResults.zmos;
        
                var dpperResults = SGP4.dpper(satrec, inclm, satrec.init, satrec.ecco, satrec.inclo, satrec.nodeo, satrec.argpo, satrec.mo, satrec.operationmode);
                
                satrec.ecco = dpperResults.ep;
                satrec.inclo = dpperResults.inclp;
                satrec.nodeo = dpperResults.nodep;
                satrec.argpo = dpperResults.argpp;
                satrec.mo = dpperResults.mp;
        
                argpm  = 0.0;
                nodem  = 0.0;
                mm     = 0.0;
                var dndt;
                
                var dsinitResults = SGP4.dsinit(whichconst, cosim, emsq, satrec.argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, satrec.t, tc, satrec.gsto, satrec.mo, satrec.mdot, satrec.no, satrec.nodeo, satrec.nodedot, xpidot, z1, z3, z11, z13, z21, z23, z31, z33, satrec.ecco, eccsq, em, argpm, inclm, mm, nm, nodem, satrec.irez, satrec.atime, satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222 , satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232, satrec.d5421, satrec.d5433, satrec.dedt, satrec.didt, satrec.dmdt, satrec.dnodt, satrec.domdt, satrec.del1, satrec.del2, satrec.del3, satrec.xfact, satrec.xlamo, satrec.xli, satrec.xni);
                
                em              = dsinitResults.em;
                argpm           = dsinitResults.argpm;
                inclm           = dsinitResults.inclm;
                mm              = dsinitResults.mm;
                nm              = dsinitResults.nm;
                nodem           = dsinitResults.nodem;
                satrec.irez     = dsinitResults.irez;
                satrec.atime    = dsinitResults.atime;
                satrec.d2201    = dsinitResults.d2201;
                satrec.d2211    = dsinitResults.d2211;
                satrec.d3210    = dsinitResults.d3210;
                satrec.d3222    = dsinitResults.d3222;
                satrec.d4410    = dsinitResults.d4410;
                satrec.d4422    = dsinitResults.d4422;
                satrec.d5220    = dsinitResults.d5220;
                satrec.d5232    = dsinitResults.d5232;
                satrec.d5421    = dsinitResults.d5421;
                satrec.d5433    = dsinitResults.d5433;
                satrec.dedt     = dsinitResults.dedt;
                satrec.didt     = dsinitResults.didt;
                satrec.dmdt     = dsinitResults.dmdt;
                dndt            = dsinitResults.dndt;
                satrec.dnodt    = dsinitResults.dnodt;
                satrec.domdt    = dsinitResults.domdt;
                satrec.del1     = dsinitResults.del1;
                satrec.del2     = dsinitResults.del2;
                satrec.del3     = dsinitResults.del3;
                satrec.xfact    = dsinitResults.xfact;
                satrec.xlamo    = dsinitResults.xlamo;
                satrec.xli      = dsinitResults.xli;
                satrec.xni      = dsinitResults.xni;
                
            }
    
            // ----------- set variables if not deep space -----------
            if (satrec.isimp !== 1) {
                cc1sq = satrec.cc1 * satrec.cc1;
                satrec.d2 = 4.0 * ao * tsi * cc1sq;
                temp = satrec.d2 * tsi * satrec.cc1 / 3.0;
                satrec.d3 = (17.0 * ao + sfour) * temp;
                satrec.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * satrec.cc1;
                satrec.t3cof = satrec.d2 + 2.0 * cc1sq;
                satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 * (12.0 * satrec.d2 + 10.0 * cc1sq));
                satrec.t5cof = 0.2 * (3.0 * satrec.d4 + 12.0 * satrec.cc1 * satrec.d3 + 6.0 * satrec.d2 * satrec.d2 + 15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq));
            }
        }
    
        SGP4.sgp4(satrec, 0.0);
        satrec.init = 'n';
    
        return true;
    },
    sgp4: function(satrec, tsince, whichconst) {
        'use strict';
        
        var am, axnl, aynl, betal, cosim, sinim, cosomm, sinomm, cnod, snod, cos2u, sin2u, coseo1, sineo1, cosi, sini, cosip, sinip, cosisq, cossu, sinsu, cosu, sinu, delm, delomg, dndt, eccm, emsq, ecose, el2, eo1, eccp, esine, argpm, argpp, omgadf, pl, r, v, rtemsq, rdotl, rl, rvdot, rvdotl, su, t2, t3, t4, tc, tem5, temp, temp1, temp2, tempa, tempe, templ, u, ux, uy, uz, vx, vy, vz, inclm, mm, nm, nodem, xinc, xincp, xl, xlm, mp, xmdf, xmx, xmy, nodedf, xnode, nodep, np;
        
        whichconst = typeof whichconst !== 'undefined' ? whichconst : SGP4.getgravconst('wgs84');
        
        var mrt = 0.0;
        if (whichconst === null) {
            whichconst = satrec.whichconst;
        }

        var temp4 =   1.5e-12;
        var twopi = 2.0 * Math.PI;
        var x2o3  = 2.0 / 3.0;

        var tumin = whichconst.tumin;
        var mu = whichconst.mu;
        var radiusearthkm = whichconst.radiusearthkm;
        var xke = whichconst.xke;
        var j2 = whichconst.j2;
        var j3 = whichconst.j3;
        var j4 = whichconst.j4;
        var j3oj2 = whichconst.j3oj2;
        
        var vkmpersec = radiusearthkm * xke/60.0;
        
        satrec.t     = tsince;
        satrec.error = 0;
        satrec.error_message = null;
        
        xmdf    = satrec.mo + satrec.mdot * satrec.t;
        var argpdf  = satrec.argpo + satrec.argpdot * satrec.t;
        nodedf  = satrec.nodeo + satrec.nodedot * satrec.t;
        argpm   = argpdf;
        mm      = xmdf;
        t2      = satrec.t * satrec.t;
        nodem   = nodedf + satrec.nodecf * t2;
        tempa   = 1.0 - satrec.cc1 * satrec.t;
        tempe   = satrec.bstar * satrec.cc4 * satrec.t;
        templ   = satrec.t2cof * t2;
        
        if (satrec.isimp !== 1) {
            delomg = satrec.omgcof * satrec.t;
            //sgp4fix use mutliply for speed instead of pow
            var delmtemp =  1.0 + satrec.eta * Math.cos(xmdf);
            delm   = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo);
            temp   = delomg + delm;
            mm     = xmdf + temp;
            argpm  = argpdf - temp;
            t3     = t2 * satrec.t;
            t4     = t3 * satrec.t;
            tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 - satrec.d4 * t4;
            tempe  = tempe + satrec.bstar * satrec.cc5 * (Math.sin(mm) - satrec.sinmao);
            templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof + satrec.t * satrec.t5cof);
        }
        
        nm    = satrec.no;
        var em    = satrec.ecco;
        inclm = satrec.inclo;
        
        if (satrec.method === 'd') {
            tc = satrec.t;
            
            var dspaceResults = SGP4.dspace(satrec.irez, satrec.d2201, satrec.d2211, satrec.d3210, satrec.d3222, satrec.d4410, satrec.d4422, satrec.d5220, satrec.d5232, satrec.d5421, satrec.d5433, satrec.dedt, satrec.del1, satrec.del2, satrec.del3, satrec.didt, satrec.dmdt, satrec.dnodt, satrec.domdt, satrec.argpo, satrec.argpdot, satrec.t, tc, satrec.gsto, satrec.xfact, satrec.xlamo, satrec.no, satrec.atime, em, argpm, inclm, satrec.xli, mm, satrec.xni, nodem, nm);
            
            var atime   = dspaceResults.atime;
            em          = dspaceResults.em;
            argpm       = dspaceResults.argpm;
            inclm       = dspaceResults.inclm;
            var xli     = dspaceResults.xli;
            mm          = dspaceResults.mm;
            var xni     = dspaceResults.xni;
            nodem       = dspaceResults.nodem;
            dndt        = dspaceResults.dndt;
            nm          = dspaceResults.nm;
        }
        
        if (nm <= 0.0) {
            satrec.error_message = 'mean motion ' + nm + ' is less than zero';
            satrec.error = 2;
            return [false, false];
        }
        
        am = Math.pow((xke / nm),x2o3) * tempa * tempa;
        nm = xke / Math.pow(am, 1.5);
        em = em - tempe;
        
        if (em >= 1.0 || em < -0.001) { // || (am < 0.95)
            satrec.error_message = 'mean eccentricity ' + em + ' not within range 0.0 <= e < 1.0';
            satrec.error = 1;
            // sgp4fix to return if there is an error in eccentricity
            return [false, false];
        }
        
        // sgp4fix fix tolerance to avoid a divide by zero
        if (em < 1.0e-6) {
            em  = 1.0e-6;
        }
        mm     = mm + satrec.no * templ;
        xlm    = mm + argpm + nodem;
        emsq   = em * em;
        temp   = 1.0 - emsq;
        
        nodem  = SGP4.fmod(nodem, twopi);
        argpm  = argpm % twopi;
        xlm    = xlm % twopi;
        mm     = (xlm - argpm - nodem) % twopi;
        
        sinim = Math.sin(inclm);
        cosim = Math.cos(inclm);
        
        var ep     = em;
        xincp  = inclm;
        argpp  = argpm;
        nodep  = nodem;
        mp     = mm;
        sinip  = sinim;
        cosip  = cosim;
        if (satrec.method === 'd') {
            var dpperResults = SGP4.dpper(satrec, satrec.inclo, 'n', ep, xincp, nodep, argpp, mp, satrec.operationmode);
            
            ep      = dpperResults.ep;
            xincp   = dpperResults.inclp;
            nodep   = dpperResults.nodep;
            argpp   = dpperResults.argpp;
            mp      = dpperResults.mp;

            if (xincp < 0.0) {
                 xincp  = -xincp;
                 nodep = nodep + Math.PI;
                 argpp  = argpp - Math.PI;
            }
            if (ep < 0.0 || ep > 1.0) {
                 satrec.error_message = 'perturbed eccentricity ' + ep + ' not within range 0.0 <= e <= 1.0'
                 satrec.error = 3;
                 // sgp4fix add return
                 return [false, false];
            }
        }
        
        if (satrec.method == 'd') {
            sinip =  Math.sin(xincp);
            cosip =  Math.cos(xincp);
            satrec.aycof = -0.5*j3oj2*sinip;
            // sgp4fix for divide by zero for xincp = 180 deg
            if (Math.abs(cosip+1.0) > 1.5e-12) {
                satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip);
            } else {
                satrec.xlcof = -0.25 * j3oj2 * sinip * (3.0 + 5.0 * cosip) / temp4;
            }
        }
        
        axnl = ep * Math.cos(argpp);
        temp = 1.0 / (am * (1.0 - ep * ep));
        aynl = ep* Math.sin(argpp) + temp * satrec.aycof;
        xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl;
        
        u    = (xl - nodep) % twopi;
        eo1  = u;
        tem5 = 9999.9;
        var ktr = 1;
        
        while (Math.abs(tem5) >= 1.0e-12 && ktr <= 10) {
            sineo1 = Math.sin(eo1);
            coseo1 = Math.cos(eo1);
            tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if (Math.abs(tem5) >= 0.95) {
                tem5 = (tem5 > 0.0 ? 0.95 : -0.95);
            }
            eo1    = eo1 + tem5;
            ktr = ktr + 1;
        }
        
        ecose = axnl*coseo1 + aynl*sineo1;
        esine = axnl*sineo1 - aynl*coseo1;
        el2   = axnl*axnl + aynl*aynl;
        pl    = am*(1.0-el2);
        
        if (pl < 0.0) {
            satrec.error_message = 'semilatus rectum ' + pl + ' is less than zero';
            satrec.error = 4;
            // sgp4fix add return
            return [false, false];
        } else {
            rl     = am * (1.0 - ecose);
            rdotl  = Math.sqrt(am) * esine/rl;
            rvdotl = Math.sqrt(pl) / rl;
            betal  = Math.sqrt(1.0 - el2);
            temp   = esine / (1.0 + betal);
            sinu   = am / rl * (sineo1 - aynl - axnl * temp);
            cosu   = am / rl * (coseo1 - axnl + aynl * temp);
            su     = Math.atan2(sinu, cosu);
            sin2u  = (cosu + cosu) * sinu;
            cos2u  = 1.0 - 2.0 * sinu * sinu;
            temp   = 1.0 / pl;
            temp1  = 0.5 * j2 * temp;
            temp2  = temp1 * temp;
    
            //  -------------- update for short period periodics ------------
            if (satrec.method == 'd') {
                cosisq = cosip * cosip;
                satrec.con41  = 3.0*cosisq - 1.0;
                satrec.x1mth2 = 1.0 - cosisq;
                satrec.x7thm1 = 7.0*cosisq - 1.0;
            }
    
            mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + 0.5 * temp1 * satrec.x1mth2 * cos2u;
            su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u;
            xnode = nodep + 1.5 * temp2 * cosip * sin2u;
            xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
            var mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / xke;
            rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u + 1.5 * satrec.con41) / xke;
    
            //  --------------------- orientation vectors -------------------
            sinsu =  Math.sin(su);
            cossu =  Math.cos(su);
            snod  =  Math.sin(xnode);
            cnod  =  Math.cos(xnode);
            sini  =  Math.sin(xinc);
            cosi  =  Math.cos(xinc);
            xmx   = -snod * cosi;
            xmy   =  cnod * cosi;
            ux    =  xmx * sinsu + cnod * cossu;
            uy    =  xmy * sinsu + snod * cossu;
            uz    =  sini * sinsu;
            vx    =  xmx * cossu - cnod * sinsu;
            vy    =  xmy * cossu - snod * sinsu;
            vz    =  sini * cossu;
    
            // --------- position and velocity (in km and km/sec) ----------
            var _mr = mrt * radiusearthkm;
            r = {x: 0.0, y: 0.0, z: 0.0};
            r["x"] = _mr * ux;
            r["y"] = _mr * uy;
            r["z"] = _mr * uz;
            v = {x: 0.0, y: 0.0, z: 0.0};
            v["x"] = (mvt * ux + rvdot * vx) * vkmpersec;
            v["y"] = (mvt * uy + rvdot * vy) * vkmpersec;
            v["z"] = (mvt * uz + rvdot * vz) * vkmpersec;
        }
        
        if (mrt < 1.0) {
            satrec.error_message = 'mrt ' + mrt + ' is less than 1.0 indicating the satellite has decayed';
            satrec.error = 6;
            return [false, false];
        }

        return {
            position: r, 
            velocity: v
        };
    },
    gstime: function(jdut1) {
        var tut1 = (jdut1 - 2451545.0) / 36525.0;
        var temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841; // sec
        temp = (temp * SGP4.deg2rad / 240.0) % SGP4.twopi; // 360/86400 = 1/240, to deg, to rad
    
        // ------------------------ check quadrants ---------------------
        if (temp < 0.0)
            temp += SGP4.twopi;
    
        return temp;
    },
    getgravconst: function(whichconst) {
        if (whichconst == 'wgs72old') {
            var mu     = 398600.79964;        //  in km3 / s2
            var radiusearthkm = 6378.135;     //  km
            var xke    = 0.0743669161;
            var tumin  = 1.0 / xke;
            var j2     =   0.001082616;
            var j3     =  -0.00000253881;
            var j4     =  -0.00000165597;
            var j3oj2  =  j3 / j2;
           //  ------------ wgs-72 constants ------------
        } else if (whichconst == 'wgs72') {
            var mu     = 398600.8;            //  in km3 / s2
            var radiusearthkm = 6378.135;     //  km
            var xke    = 60.0 / Math.sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
            var tumin  = 1.0 / xke;
            var j2     =   0.001082616;
            var j3     =  -0.00000253881;
            var j4     =  -0.00000165597;
            var j3oj2  =  j3 / j2;
        } else if (whichconst == 'wgs84') {
            // ------------ wgs-84 constants ------------
            var mu     = 398600.5;            //  in km3 / s2
            var radiusearthkm = 6378.137;     //  km
            var xke    = 60.0 / Math.sqrt(radiusearthkm*radiusearthkm*radiusearthkm/mu);
            var tumin  = 1.0 / xke;
            var j2     =   0.00108262998905;
            var j3     =  -0.00000253215306;
            var j4     =  -0.00000161098761;
            var j3oj2  =  j3 / j2;
        }
        return {
            tumin: tumin,
            mu: mu,
            radiusearthkm: radiusearthkm,
            xke: xke,
            j2: j2,
            j3: j3,
            j4: j4,
            j3oj2: j3oj2,
        };
    },
    wgs72: function() {
        return SGP4.getgravconst('wgs72');
    },
    wgs84: function() {
        return SGP4.getgravconst('wgs84');
    },
    twoline2rv: function(line1, line2, gravconst) {
        'use strict';
        
        var opsmode = 'i';
        
        var deg2rad  =  Math.PI / 180.0; //0.0174532925199433
        var xpdotp   =  1440.0 / (2.0 * Math.PI); //229.1831180523293
        
        var tumin = gravconst.tumin;
        
        var satrec = {};
        satrec.error = 0;
        satrec.whichconst = gravconst;
        
        //Line 1
        satrec.satnum = line1.substring(2, 7);
        satrec.epochyr = parseInt(line1.substring(18, 20), 10);
        satrec.epochdays = parseFloat(line1.substring(20, 32));
        satrec.ndot = parseFloat(line1.substring(33, 43));
        satrec.nddot = parseFloat("." + parseInt(line1.substring(44, 50), 10) + "E" + line1.substring(50, 52));
        satrec.bstar = parseFloat(line1.substring(53, 54) + '.' + line1.substring(54, 59)); //parseFloat("." + parseInt(line1.substring(53, 59), 10) + "E" + line1.substring(59, 61));
        var ibexp = parseInt(line1.substring(59, 61), 10);
        satrec.bstar = satrec.bstar * Math.pow(10, ibexp);
        
        //Line 2
        satrec.inclo = parseFloat(line2.substring(8, 16));
        satrec.nodeo = parseFloat(line2.substring(17, 25));
        satrec.ecco = parseFloat("." + line2.substring(26, 33));
        satrec.argpo = parseFloat(line2.substring(34, 42));
        satrec.mo = parseFloat(line2.substring(43, 51));
        satrec.no = parseFloat(line2.substring(52, 63));
        
        satrec.no = satrec.no / xpdotp; //   rad/min
        satrec.ndot = satrec.ndot  / (xpdotp * 1440.0);  // ? * minperday
        satrec.nddot= satrec.nddot / (xpdotp * 1440.0 * 1440);
        
        satrec.inclo = satrec.inclo  * deg2rad;
        satrec.nodeo = satrec.nodeo  * deg2rad;
        satrec.argpo = satrec.argpo  * deg2rad;
        satrec.mo    = satrec.mo     * deg2rad;
        
        satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
        satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;
        
        var year = 0;
        if (satrec.epochyr < 57) {
            year = satrec.epochyr + 2000;
        } else {
            year = satrec.epochyr + 1900;
        }
        
        var days2mdhmsResult = SGP4.days2mdhms(year, satrec.epochdays);
        var mon, day, hr, minute, sec;
        mon = days2mdhmsResult.mon;
        day = days2mdhmsResult.day;
        hr = days2mdhmsResult.hr;
        minute = days2mdhmsResult.minute;
        sec = days2mdhmsResult.sec;
        
        satrec.jdsatepoch = SGP4.jday(year, mon, day, hr, minute, sec);
        
        SGP4.sgp4init(gravconst, opsmode, satrec.satnum, satrec.jdsatepoch-2433281.5, satrec.bstar, satrec.ecco, satrec.argpo, satrec.inclo, satrec.mo, satrec.no, satrec.nodeo, satrec);
        
        return satrec;
    },
    days2mdhms: function(year, days) {
        var lmonth = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
    
        var dayofyr = Math.floor(days);
        //  ----------------- find month and day of month ----------------
        if (year % 4 == 0)
            lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
    
        var i = 1;
        var inttemp = 0;
        while (dayofyr > inttemp + lmonth[i-1] && i < 12) {
            inttemp = inttemp + lmonth[i-1];
            i += 1;
        }
    
        var mon = i;
        var day = dayofyr - inttemp;
    
        //  ----------------- find hours minutes and seconds -------------
        var temp = (days - dayofyr) * 24.0;
        var hr   = Math.floor(temp);
        temp = (temp - hr) * 60.0;
        var minute  = Math.floor(temp);
        var sec  = (temp - minute) * 60.0;
    
        return {
            mon: mon,
            day: day,
            hr: hr,
            minute: minute,
            sec: sec
        };
    },
    jday: function(year, mon, day, hr, minute, sec) {
        return (367.0 * year - Math.floor((7 * (year + Math.floor((mon + 9) / 12.0))) * 0.25) + Math.floor( 275 * mon / 9.0 ) + day + 1721013.5 + ((sec / 60.0 + minute) / 60.0 + hr) / 24.0);
    },
    propogate: function(satrec, year, month, day, hour, minute, second) {
        var j = SGP4.jday(year, month, day, hour, minute, second);
        var m = (j - satrec.jdsatepoch) * 1440;
        return SGP4.sgp4(satrec, m);
    },
    // All functions below are taken from satellite-js found at (https://github.com/shashwatak/satellite-js), created by Shashwat Kandadai and released under the MIT license.
    // Some modifications have been made to the original code.
    /*
     * satellite-js v1.1
     * (c) 2013 Shashwat Kandadai and UCSC
     * https://github.com/shashwatak/satellite-js
     * License: MIT
     */
    eciToGeodetic: function(eci_coords, gmst) {
        'use strict';
        // http://www.celestrak.com/columns/v02n03/
        var a   = 6378.137;
        var b   = 6356.7523142;
        var R   = Math.sqrt( (eci_coords["x"]*eci_coords["x"]) + (eci_coords["y"]*eci_coords["y"]) );
        var f   = (a - b)/a;
        var e2  = ((2*f) - (f*f));
        var longitude = Math.atan2(eci_coords["y"], eci_coords["x"]) - gmst;
        var kmax = 20;
        var k = 0;
        var latitude = Math.atan2(eci_coords["z"],
                       Math.sqrt(eci_coords["x"]*eci_coords["x"] +
                                    eci_coords["y"]*eci_coords["y"]));
        var C;
        while (k < kmax){
            C = 1 / Math.sqrt( 1 - e2*(Math.sin(latitude)*Math.sin(latitude)) );
            latitude = Math.atan2 (eci_coords["z"] + (a*C*e2*Math.sin(latitude)), R);
            k += 1;
        }
        var height = (R/Math.cos(latitude)) - (a*C);
        
        var velocity = Math.sqrt(398600.8 / (height + 6378.135)); // Velocity in kilemeters per second. Multiply by 3600 for kilometers per hour.
        
        return { longitude : longitude, latitude : latitude, height : height, velocity : velocity };
    },
    eciToEcf: function(eci_coords, gmst){
        'use strict';
        // ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
        //
        // [X]     [C -S  0][X]
        // [Y]  =  [S  C  0][Y]
        // [Z]eci  [0  0  1][Z]ecf
        //
        //
        // Inverse:
        // [X]     [C  S  0][X]
        // [Y]  =  [-S C  0][Y]
        // [Z]ecf  [0  0  1][Z]eci
    
        var X = (eci_coords["x"] * Math.cos(gmst))    + (eci_coords["y"] * Math.sin(gmst));
        var Y = (eci_coords["x"] * (-Math.sin(gmst))) + (eci_coords["y"] * Math.cos(gmst));
        var Z =  eci_coords["z"];
        return { x : X, y : Y, z : Z };
    },
    ecfToEci: function(ecf_coords, gmst){
        'use strict';
        // ccar.colorado.edu/ASEN5070/handouts/coordsys.doc
        //
        // [X]     [C -S  0][X]
        // [Y]  =  [S  C  0][Y]
        // [Z]eci  [0  0  1][Z]ecf
        //
        var X = (ecf_coords["x"] * Math.cos(gmst))    - (ecf_coords["y"] * Math.sin(gmst));
        var Y = (ecf_coords["x"] * (Math.sin(gmst)))  + (ecf_coords["y"] * Math.cos(gmst));
        var Z =  ecf_coords["z"];
        return { x : X, y : Y, z : Z };
    },
    geodeticToEcf: function(geodetic_coords){
        'use strict';
        var longitude   = geodetic_coords["longitude"];
        var latitude    = geodetic_coords["latitude"];
        var height      = geodetic_coords["height"];
        var a           = 6378.137;
        var b           = 6356.7523142;
        var f           = (a - b)/a;
        var e2          = ((2*f) - (f*f));
        var normal      = a / Math.sqrt( 1 - (e2*(Math.sin(latitude)*Math.sin(latitude))));
    
        var X           = (normal + height) * Math.cos (latitude) * Math.cos (longitude);
        var Y           = (normal + height) * Math.cos (latitude) * Math.sin (longitude);
        var Z           = ((normal*(1-e2)) + height) * Math.sin (latitude);
        return { x : X, y : Y, z : Z };
    },
    topocentric: function(observer_coords, satellite_coords){
        // http://www.celestrak.com/columns/v02n02/
        // TS Kelso's method, except I'm using ECF frame
        // and he uses ECI.
        //
        'use strict';
        var longitude   = observer_coords["longitude"];
        var latitude    = observer_coords["latitude"];
        var height      = observer_coords["height"];
    
        var observer_ecf = SGP4.geodeticToEcf(observer_coords);
    
        var rx      = satellite_coords["x"] - observer_ecf["x"];
        var ry      = satellite_coords["y"] - observer_ecf["y"];
        var rz      = satellite_coords["z"] - observer_ecf["z"];
    
        var top_s   = ( (Math.sin(latitude) * Math.cos(longitude) * rx) +
                      (Math.sin(latitude) * Math.sin(longitude) * ry) -
                      (Math.cos(latitude) * rz));
        var top_e   = ( -Math.sin(longitude) * rx) + (Math.cos(longitude) * ry);
        var top_z   = ( (Math.cos(latitude)*Math.cos(longitude)*rx) +
                      (Math.cos(latitude)*Math.sin(longitude)*ry) +
                      (Math.sin(latitude)*rz));
        return { top_s : top_s, top_e : top_e, top_z : top_z };
    },
    topocentricToLookAngles: function(topocentric){
        'use strict';
        var top_s = topocentric["top_s"];
        var top_e = topocentric["top_e"];
        var top_z = topocentric["top_z"];
        var range_sat    = Math.sqrt((top_s*top_s) + (top_e*top_e) + (top_z*top_z));
        var El      = Math.asin (top_z/range_sat);
        var Az      = Math.atan2 (-top_e, top_s) + Math.PI;
        return { azimuth : Az, elevation : El, range_sat : range_sat };
    },
    degreesLong: function(radians){
        'use strict';
        var degrees = (radians/Math.PI*180) % (360);
        if (degrees > 180){
            degrees = 360 - degrees;
        }
        else if (degrees < -180){
            degrees = 360 + degrees;
        }
        return degrees;
    },
    degreesLat: function(radians){
        'use strict';
        if (radians > Math.PI/2 || radians < (-Math.PI/2)){
            return "Err";
        }
        var degrees = (radians/Math.PI*180);
        if (degrees < 0){
            degrees = degrees;
        }
        else{
            degrees = degrees;
        }
        return degrees;
    },
    gstimeFromDate: function(year, mon, day, hr, minute, sec) {
        var jDay = SGP4.jday(year, mon, day, hr, minute, sec);
        return SGP4.gstime (jDay);
    },
    gstime: function(jd) {
        var tut1 = (jd - 2451545.0) / 36525.0;
        var temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 +
                (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841;  //#  sec
        temp = (temp * SGP4.deg2rad / 240.0) % SGP4.twopi; // 360/86400 = 1/240, to deg, to rad
    
        //  ------------------------ check quadrants ---------------------
        if (temp < 0.0){
            temp += SGP4.twopi;
        }
        return temp;
    },
};
module.exports = SGP4;
