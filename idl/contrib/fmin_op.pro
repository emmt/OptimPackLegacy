; $Id: fmin_op.pro,v 1.16 2007/04/02 15:12:12 mugnier Exp $
; Copyright (c) L. Mugnier, ONERA, 2007.
;
; This software is copyright ONERA.
; The author is Laurent Mugnier.
; E-mail: mugnier at onera.fr
;
; This software is a computer program whose purpose is to be an IDL frontend
; to the OptimPack package by Eric Thiébaut (CRAL) that is (almost) a
; plug-and-play replacement for previous minimization engines used at ONERA.
;
; This software is governed by the CeCILL-C license under French law and
; abiding by the rules of distribution of free software. You can use, modify
; and/ or redistribute the software under the terms of the CeCILL-C license as
; circulated by CEA, CNRS and INRIA at the following URL
; "http://www.cecill.info".
;
; As a counterpart to the access to the source code and  rights to copy,
; modify and redistribute granted by the license, users are provided only
; with a limited warranty  and the software's author,  the holder of the
; economic rights,  and the successive licensors  have only  limited
; liability. 
;
; In this respect, the user's attention is drawn to the risks associated
; with loading,  using,  modifying and/or developing or reproducing the
; software by the user in light of its specific status of free software,
; that may mean  that it is complicated to manipulate,  and  that  also
; therefore means  that it is reserved for developers  and  experienced
; professionals having in-depth computer knowledge. Users are therefore
; encouraged to load and test the software's suitability as regards their
; requirements in conditions enabling the security of their systems and/or 
; data to be ensured and,  more generally, to use and operate it in the 
; same conditions as regards security. 
;
; The fact that you are presently reading this means that you have had
; knowledge of the CeCILL-C license and that you accept its terms.
;
PRO FMIN_OP, x, fx, FUNC = func, INIT = init, $
             WS = ws, GRADIENT = gradient, $ ; entrees/sorties necessaires
             ERRF = errf, $     ; sortie utile
             JOB = job, $       ; sortie optionnelle
             ITER = iter, $     ; sortie optionnelle
             REINIT = reinit, $ ; entree optionnelle
             CONV_THRESHOLD = conv_threshold, $ ; entree optionnelle
             ITMAX = itmax, $   ; entree optionnelle
             _REF_EXTRA = extra, VISU = visu, $ ; entrees optionnelles
             ACTIVE_SET = active_set, $ ; entree pour contrainte de positivite
             XMIN = xmin, XMAX = xmax, $
             MEMORY = memory, $ ; entree optionnelle
             LIBRARY = library, ALLTHEWAY = alltheway, $
             VERSION = version, HELP = help
;+
;NOM :
;   FMIN_OP - Minimisation de critère par appel à OptimPack (VMLM avec Bornes)
;   
;CATEGORIE :
;   Mathematics Routines
;
;SYNTAXE :
;   PRO FMIN_OP, x [, fx] , FUNC = func, INIT = init, $
;                WS = ws, GRADIENT = gradient, $;entrees/sorties necessaires
;                [ERRF = errf, $]                    ; sortie utile
;                [JOB = job, $]                      ; sortie optionnelle
;                [ITER = iter, $]                    ; sortie optionnelle
;                [REINIT = reinit, $]                ; entree optionnelle
;                [CONV_THRESHOLD = conv_threshold, $]; entree optionnelle
;                [ITMAX = itmax, $]                  ; entree optionnelle
;                [_REF_EXTRA = extra][, VISU = visu, $];entrees optionnelles
;                [ACTIVE_SET = active_set, $]        ; entree optionnelle
;                [XMIN = xmin, ][XMAX = xmax, $]     ; entrees optionnelles
;                [MEMORY = memory, $]                ; entree optionnelle
;                [LIBRARY=library, $]                ; entree optionnelle
;                [/VERSION] [, /HELP]
;
;DESCRIPTION :
;   FMIN_OP minimise une fonction de N variables par une méthode de type VMLM
;   avec Bornes. C'est un ``wrapper'' ou frontal des routines op_init,
;   op_vmlmb_setup et op_vmlmb développées par Eric Thiébaut (package
;   OptimPack) destiné à remplacer FMIN_GPA (syntaxe identique). FMIN_OP
;   doit être appelé dans une boucle (car effectue une itération de
;   minimisation) sauf avec le mot-cle /ALLTHEWAY.
;   
;   ARGUMENTS :
;
;     x             : (entrée/sortie) valeur courante de l'estimée du point
;                     minimum (en 1ère entrée, point de départ de la
;                     minimisation).
;
;    fx             : (entrée/sortie optionnelle) nom de variable, qui
;                     contient en sortie la valeur du critère au point x.
;                     Utilisable en entrée comme valeur précédente de f(x)
;                     pour calcul de ERRF.
;
;    FUNC           : (entrée) nom de la fonction qui calcule le critère et
;                     son gradient.
;
;    INIT           : (entrée) à mettre a 1 lors du premier appel puis à 0
;                     pour toute la suite de la minimisation. 
;                     Inutile avec le mot-clé /ALLTHEWAY.
;
;    WS             : (entrée/sortie) nom de variable qui reçoit en sortie et
;                     pour les itérations suivantes le ``WorkSpace'' pour
;                     VMLM-B (cf. op_vmlmb.pro). 
;                     Inutile avec le mot-clé /ALLTHEWAY. 
;
;    GRADIENT       : (entrée/sortie) nom de variable qui reçoit en sortie et
;                     pour les itérations suivantes le gradient courant.
;                     Inutile avec le mot-clé /ALLTHEWAY.
;
;    ERRF           : (sortie optionnelle) utile contient la valeur du test de
;                     convergence sur FUNC(x), qui est l'évolution du critère
;                     FUNC(x) entre les 2 dernières itérations, relativement à
;                     la valeur moyenne du critère.
;
;
;    JOB            : (sortie optionnelle) utile si on veut l'avis de op_vmlmb
;                     sur la convergence. Contient la valeur de JOB rendue
;                     par op_vmlmb (voir ce programme). Doit valoir 2 si
;                     itération bien terminée, 3 si op_vmlmb considère qu'il a
;                     convergé (cf CONV_THRESHOLD), et 4 s'il n'arrive plus à
;                     avancer.
;
;    ITER           : (sortie optionnelle) contient le nombre d'itérations
;                     effectuées. 
;
;    REINIT         : (entrée optionnelle) permet de réinitialiser
;                     l'algorithme depuis le point courant si présent et
;                     non nul (méthode pour 1ère descente = gradient simple). 
;                     Par défaut on ne réinitialise jamais.
;
;    CONV_THRESHOLD : (entrée optionnelle) utile si on veut l'avis de op_vmlmb
;                     sur la convergence. Si présent, utilisé par op_vmlmb
;                     comme seuil sur l'évolution du critère en relatif.
;                     Mot-clé indispensable avec /ALLTHEWAY.
;
;    ITMAX          : (entrée optionnelle) nombre maximum d'itérations. A
;                     n'utiliser que si l'on ne veut pas minimiser
;                     complètement le critère. 
;
;    _REF_EXTRA     : (entrée) permet de passer des variables et mots-clés à
;                     FUNC. 
;
;    VISU           : (entrée optionnelle) visu d'infos pendant la
;                     minimisation, toutes les VISU itérations (ainsi qu'à
;                     l'itération 0 et si job GT 2 càd convergence atteinte).
;                     Les infos visualisées sont job, le nb d'itérations
;                     iter_op, le nb d'évaluation du critère et de son
;                     gradient eval_op, la valeur du critère f(x) et la valeur
;                     du test de convergence (cf mot-clé CONV_THRESHOLD).
;    
;    ACTIVE_SET     : (entrée/sortie optionnelle) pour demander une
;                     minimisation sous contrainte de positivité ; tableau de
;                     type BYTE de la taille de x qui doit valoir 1B lors du
;                     premier appel, et vaut 1B sur les points actifs ensuite.
;
;    MEMORY         : (entrée optionnelle) détermine la taille mémoire
;                     utilisée pour approximer l'inverse du Hessien par la
;                     méthode "Variable metric with Limited Memory" i.e. BFGS
;                     à mémoire limitée. MEMORY=5 par défaut (et est
;                     proportionnelle au nombre de gradients conservés en
;                     mémoire). 
; 
;   LIBRARY=library : (entrée optionnelle) nom complet de la bibliothèque
;                     OptimPack_IDL${OSTYPE} (avec répertoire mais sans le .so,
;                     comme requis par OP_INIT). Exemple = OptimPack_IDLlinux
;                     sous linux, OptimPack_IDLsolaris sous solaris.
;                     Sous Unix ce mot-clé est inutile car cette bibliothèque
;                     est trouvée automatiquement dès qu'elle est dans le
;                     !PATH (cf. whereis.pro). Sur les machines 64 bits,
;                     c'est OP_INIT qui ajoutera _64 au nom de librairie
;                     à charger, donnant OptimPack_IDL${OSTYPE}_64.so (cf doc
;                     de OP_INIT). 
;
;   /ALLTHEWAY      : (entrée optionnelle) si ce mot-clé est présent et
;                     non nul, FMIN_OP itère jusqu'à ce que (JOB GE 3)
;                     i.e. fin de convergence au lieu de (JOB GE 2) par défaut
;                     i.e. fin d'une itération.
;                     Très utile si on ne veut rien afficher entre 2 itérations.
;                     Avec /ALLTHEWAY, les mots-clés INIT, WS et GRADIENT ne
;                     sont plus obligatoires.
;                     CONV_THRESHOLD devient indispensable.
;   
;   /VERSION        : (entrée) affichage de la version avant l'exécution.
;                     Si (VERSION GE 2) alors VERSION est également passé à
;                     FUNC pour avoir sa version (ce qui suppose que FUNC
;                     accepte le mot-clé VERSION).
;   
;   /HELP           : (entrée) affichage de la syntaxe et sortie du programme.
;
;DIAGNOSTIC D'ERREUR :
;   
;EXEMPLE :
;
;   active_set = byte(x) * 0B + 1B ; où x est le "guess" initial
;   criterion_value = -1.
;   FMIN_OP, x, criterion_value, FUNC = 'criterion_computation_function', $
;            /ALLTHEWAY, CONV_THRESHOLD = (machar()).eps, $
;            ITMAX = 1000L,  ACTIVE_SET = active_set, $
;
;
;VOIR AUSSI :
;   OP_INIT, OP_VMLMB_SETUP, OP_VMLMB et OP_VMLMB_MSG (Eric Thiébaut),
;   WHEREIS.
;
;AUTEUR :
;   $Author: mugnier $
;
;HISTORIQUE :
;   $Log: fmin_op.pro,v $
;   Revision 1.16  2007/04/02 15:12:12  mugnier
;   Amélioration de la doc et exemple.
;
;   Revision 1.15  2007/01/31 15:04:55  mugnier
;   Doc sur infos affichées.
;
;   Revision 1.14  2007/01/29 15:55:37  mugnier
;   VERSION peut désormais être passé à FUNC (voir doc).
;
;   Revision 1.13  2007/01/26 14:29:39  mugnier
;   Mot-clé VERSION passé à WHEREIS.
;
;   Revision 1.12  2007/01/26 11:46:13  mugnier
;   Added copyright and CeCILL-C license in source.
;
;   Revision 1.11  2007/01/23 17:43:11  mugnier
;   Correction d'un bug si le chemin de la bibliothèque OptimPack contient "../"
;
;   Revision 1.10  2007/01/19 17:50:15  mugnier
;   Correction petit bug affichage info (si on n'en veut pas).
;
;   Revision 1.9  2007/01/16 17:34:51  mugnier
;   INIT=1 n'est plus nécessaire lors d'un appel avec le mot-clé /ALLTHEWAY.
;
;   Revision 1.8  2006/10/31 15:12:10  mugnier
;   Nouveau mot-clé ITMAX.
;
;   Revision 1.7  2006/04/27 14:55:14  mugnier
;   Mot-clé /ALLTHEWAY pour minimisation complète au lieu d'une seule itération.
;
;   Revision 1.6  2006/04/27 13:12:08  mugnier
;   Doc sur library mise à jour conjointement à op_init.
;
;   Revision 1.5  2005/06/08 12:56:19  mugnier
;   La bibliothèque recherchée par défaut est désormais OptimPack_IDL${OSTYPE},
;   soit par exemple OptimPack_IDLlinux ou OptimPack_IDLsolaris.
;   Ceci permet de garder au même endroit des bibliothèques pour différents OSs.
;
;   Revision 1.4  2004/02/12 09:36:31  mugnier
;   Nouveau mot-clé REINIT pour éventuellement réinitialiser la direction
;   de descente régulièrement si stagnation constatée. *Devrait* être inutile.
;   Implémentation OK d'après E. Thiébaut mais à tester plus complètement.
;
;   Revision 1.3  2003/10/02 15:17:20  mugnier
;   - FX peut etre désormais utilisé en entrée (cf. doc) ;
;   - mots-clés XMIN et XMAX (bornes autres que 0 et infinité) ;
;   - correction d'un effet de bord sur le mot-clé VISU.
;
;   Revision 1.2  2003/04/22 15:51:46  mugnier
;   Message d'erreur si version d'IDL < 5.4 : OptimPack teste la valeur de
;     !version.memory_bits.
;   Corrigé bug si conv_threshold était un nom de variable sans variable.
;   Le 2ème argument, fx, est désormais optionnel comme avec FMIN_GPA.
;
;   Revision 1.1  2003/04/16 16:49:42  mugnier
;   Initial revision
;
;-

on_error,2
IF keyword_set(version) THEN BEGIN 
   version_inside = version
   printf, -2, '% ' + routine_courante() + $
            ': $Revision: 1.16 $, $Date: 2007/04/02 15:12:12 $'
ENDIF ELSE $
    version_inside = 0L

IF (NOT ((n_params() GT 0L) AND $
         ((arg_present(ws) AND arg_present(gradient)) $
          OR keyword_set(alltheway)))) $
    OR keyword_set(help) THEN BEGIN 
    message, 'Help wanted or incorrect syntax. Documentation :', /INFO
    doc_library, 'fmin_op'
    retall
ENDIF

IF (keyword_set(active_set)) THEN BEGIN
    projection = 1B
    IF ((size(active_set, /type) NE 1L) OR $ ; type 1 = byte
        (n_elements(active_set) NE n_elements(x))) THEN message, $
            'si positivité par projection, active_set doit être bytarr ' + $
            'et de taille de x'
    IF (n_elements(xmax) EQ 0L) THEN xmax = (machar(/double)).xmax
    IF (n_elements(xmin) EQ 0L) THEN xmin = 0.0D
ENDIF ELSE projection = 0B

IF keyword_set(visu) THEN $ ; Ne pas modifier `visu' si passé et 0
                visu_long = visu $
ELSE $
                visu_long = 2147483647L ;=2L^31-1: pas de visu par defaut. 

IF keyword_set(alltheway) THEN job_endvalue = 3L ELSE job_endvalue = 2L
IF n_elements(itmax) EQ 0L THEN itmax = 2147483647L ; 2^31-1

COMMON fmin_op, job_op, eval_op, iter_op, fx_op, errf_op, negsetcount_op

IF keyword_set(init) OR keyword_set(alltheway) THEN BEGIN ;initialisation 
    job_op = 1L
    eval_op = 0L
    iter_op = 0L
    IF (n_elements(fx) NE 0L) THEN fx_op = double(fx) ELSE fx_op = 0.0D
    errf_op = 0.0D
    
    IF NOT keyword_set(library) THEN BEGIN
        defaultlibname = 'OptimPack_IDL'+getenv('OSTYPE')+'.so'
        whereis, defaultlibname, /NOSUFFIX, chemin, /QUIET, $
                 VERSION = version_inside
        ; NB: whereis only work with Unix-like systems
        IF chemin[0] EQ '' THEN $
            message, 'bibliothèque '+defaultlibname+' introuvable !' $
        ELSE $ ;chemin sans .so :
;            library = (strsplit(chemin[0], '.', /extract))[0]
; la ligne ci-dessus deconne si chemin[0] contient "../" donc on fait ceci :                     
             library = (strsplit(chemin[0], '\.so', /extract,/regex))[0]
    ENDIF
    IF (float(!version.release) LT 5.4) THEN $
        message, 'op_init (donc FMIN_OP) ne marche pas avant IDL 5.4.'
;    print, "BIBLIOTHEQUE :", library
    op_init, library
    
    IF keyword_set(memory) THEN memory_long = long(memory) $
        ELSE memory_long = 5L
    IF (n_elements(conv_threshold) EQ 0L) THEN conv_threshold = 0.0D
    ws = op_vmlmb_setup(n_elements(x), memory_long, $
                        frtol=double(conv_threshold), fatol=0.0D)
ENDIF

IF keyword_set(reinit) THEN reinit_NOT_done = 1L ELSE  reinit_NOT_done = 0L

last_fx  = fx_op ;  fx_op=valide en entree si job GE 2

REPEAT BEGIN
    IF (job_op EQ 1L) AND keyword_set(reinit) AND reinit_NOT_done THEN BEGIN
        IF keyword_set(memory) THEN memory_long = long(memory) $
            ELSE memory_long = 5L
        IF (n_elements(conv_threshold) EQ 0L) THEN conv_threshold = 0.0D
        ws = op_vmlmb_setup(n_elements(x), memory_long, $
                            frtol=double(conv_threshold), fatol=0.0D)
        reinit_NOT_done = 0L
    ENDIF

    IF (job_op EQ 1L) THEN BEGIN ;évaluation du critère et du gradient en x
        ; projection de x à faire avant :
        IF projection THEN BEGIN
            IF ((iter_op MOD visu_long) EQ 0L) THEN $
                negative_set = where((x LT xmin) OR (x GT xmax), negsetcount_op)
            x = (temporary(x) >  xmin) < xmax
        ENDIF
        
;        last_fx = fx_op ; fx_op=dernier fx evalué qd on demande justement fx
        
        IF (version_inside GE 2) THEN BEGIN
           IF keyword_set(extra) THEN $;on veut toujours gradient=>tjrs passé
               fx_op = CALL_FUNCTION(func, x, gradient, _EXTRA=extra, $
                                     VERSION = version_inside) $
           ELSE $
               fx_op = CALL_FUNCTION(func, x, gradient, VERSION=version_inside)
        ENDIF ELSE BEGIN 
           IF keyword_set(extra) THEN $;on veut toujours gradient=>tjrs passé
               fx_op = CALL_FUNCTION(func, x, gradient, _EXTRA=extra) $ ;_STRICT
           ELSE $
               fx_op = CALL_FUNCTION(func, x, gradient)
        ENDELSE

        eval_op = eval_op + 1L
;        errf_op = 2.*abs(fx_op - last_fx)/(abs(fx_op)+abs(last_fx));valide en sortie si job=2|3
    ; f(x) est valide en sortie ssi job=2 ou 3 (cf op_vmlm.pro).
    ; donc errf est valide en sortie ssi job=2 ou 3 
    ; en effet si job=2 ou 3 en sortie de op_vmlmb on avait job=1 en entrée.
    ENDIF

    ;; Call optimizer.
    ;; NB : gradient doit etre present meme si pas actuel (ie meme si job NE 1)
    
    IF projection THEN BEGIN 
        IF ((job_op GT 1L) OR (eval_op EQ 1L)) THEN BEGIN
            ;; mise à jour de l'active set (job > 1 ou 1ere eval critere)
;            active_set = ((x GT 0.) OR (gradient LT 0.))
            active_set = ((x GT xmin) OR (gradient LT 0.0D)) AND $
                         ((x LT xmax) OR (gradient GT 0.0D))
;; On est sûr qu'on n'a pas ici de valeurs de x <0 parce que si x est modifié
;; par op_vmlmb alors forcément job=1 en sortie donc on applique les bornes
;; juste après.      
            
            IF (((iter_op+1L) MOD visu_long) EQ 0L) THEN BEGIN 
            zero_set = where((x EQ xmin) OR (x EQ xmax), count1)
            zero_and_positivegrad_set = where(active_set EQ 0B, count) 
            print, 'Nb pts < 0 : ', nbr2str(negsetcount_op), $
               '. Après proj, nb pts vérifiant Kuhn-Tucker/nuls : ', $
               nbr2str(count)+'/'+nbr2str(count1)+'.'
            ENDIF 
        ENDIF
        job_op = OP_VMLMB(x, fx_op, gradient, ws, ACTIVE = active_set)
    ENDIF ELSE $
        job_op = OP_VMLMB(x, fx_op, gradient, ws)
    
    IF (job_op GE 2L) THEN iter_op = iter_op + 1L

 ENDREP UNTIL ((job_op GE job_endvalue) OR (iter_op GT itmax))
; sortir avec job=2 ou 3, i.e. à la fin d'une iter ou à convergence,
; ou apres itmax iterations (et job GE 2).

errf_op = (2.0D*abs(last_fx - fx_op))/(abs(fx_op)+abs(last_fx))

IF (((iter_op MOD visu_long) EQ 0L) OR job_op GT 2L) THEN $
    printf,-2, 'FMIN_OP: Job=', nbr2str(job_op), $
            '; iter_op=', nbr2str(iter_op),$
            '; eval_op=', nbr2str(eval_op), '; f(x)=', $
            nbr2str(fx_op, format = '(G16.10)'), '; conv=', errf_op

;; passage des valeurs en sortie (okazou on les veut)
fx = fx_op
errf = errf_op
job = job_op
iter = iter_op

END
