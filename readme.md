
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
![GitHub all releases](https://img.shields.io/github/downloads/XinyueMa-neuro/PVIN-model-MaEtAl2023/total)
![GitHub followers](https://img.shields.io/github/followers/XinyueMa-neuro?style=social)
![GitHub watchers](https://img.shields.io/github/watchers/XinyueMa-neuro/PVIN-model-MaEtAl2023?style=social)
![GitHub Repo stars](https://img.shields.io/github/stars/XinyueMa-neuro/PVIN-model-MaEtAl2023?style=social)
![Twitter Follow](https://img.shields.io/twitter/follow/XinyueMa_neuro?style=social)


In this repository is the code associated with the following paper:

>**Ma, X.**, Miraucourt, L. S., Qiu, H., Sharif-Naeini, R., & Khadra, A. (2023). Modulation of SK channels via calcium buffering tunes intrinsic excitability of parvalbumin interneurons in neuropathic pain: A computational and experimental investigation. Journal of Neuroscience, 43(31), 5608-5622.

The code can be used to produce the time series simulations and bifurcation analysis of PVIN model.
  
---
<p>
Author: Xinyue Ma <br>
Email: xinyue.ma@mail.mcgill.ca <br>
Integrated Program in Neuroscience <br>
McGill University <br>
Montreal, QC, H3A 1A1 <br> 
Canada <br>
</p>


## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

To reference this code, please cite the journal article mentioned above.

---

<h2>I. Description</h2>
<p>
The PVIN model is a Hodgkin-Huxley type model of parvalbumin-expressing interneurons (PVINs) in the dorsal horn of the spinal cord. It consists ion channels including fast-transient sodium channel (Nav), slow delayed-rectifier potassium channel (Kv1), fast delayed-rectifier potassium channel (Kv3), high voltage-gated calcium channel (Cav), small conductance calcium-activated potassium channel (SK), and a leak current. It also uses the flux-balanced equation describing the intracellular calcium dynamics. The PVIN model is simulated using a stable step current protocol and a randomized synaptic current protocol. The code of PVIN model written in Matlab and XPP are both provided.
</p>
<p>
The model is adopted from Bischop et al. 2012 and reparametrized to fit the electrical activity of spinal dorsal horn PVINs from naive mice. It reproduces the firing pattern change of PVINs from tonic to transient following nerve injury, which is realized by reducing cytosolic calcium buffer concentration. The bifurcation analysis of PVIN model further explains how the firing pattern transits as the injection current increases, in a manner similar to that is seen in our transient firing PVIN recordings.
</p>
<p>
The code also includes a “in vivo-like” neural circuit model of A&#946 fiber-mediate nociceptive neural circuit. The circuit model is stimulated by Poisson-distributed excitatory synaptic currents representing the presynaptic inputs from the A&#946 fibers. It includes the PVIN model above, as well as another HH type model describing the excitability of a PVIN post-synaptic target: the excitatory interneuron expressing protein kinase C gamma (PKC&#947IN). The PKC&#947IN model is adapted from Medlock et al. 2022. The A&#946 fiber-like presynaptic current was applied on both the inhibitory PVIN model and the excitatory PKCγIN model, the latter of which also received inhibitory synaptic input from the PVIN model. 
</p>
<p>
<b>Reference:</b> <br>
<p>&emsp;Bischop, D. P., Orduz, D., Lambot, L., Schiffmann, S. N., & Gall, D. (2012). Control of neuronal excitability by calcium binding proteins: a new mathematical model for striatal fast-spiking interneurons. Frontiers in molecular neuroscience, 5, 78.</p>
<p>&emsp;Medlock L, Sekiguchi K, Hong S, Dura-Bernal S, Lytton WW, Prescott SA (2022) Multiscale Computer Model of the Spinal Dorsal Horn Reveals Changes in Network Processing Associated with Chronic Pain. J Neurosci 42:3133–3149.
</p>
</p>

<h2>II. Simulators</h2> 
    MATLAB Version: 9.12.0 (R2022a) <br>
    XPPAUT Version: 8.0 <br>
    Time-stamp: 2023-March-05 <br>

  
<h2>III. File list</h2>
 
### XPP
<p>
    <table style="width:100%">
        <tr>
            <th><code>PVIN_2023.ode</code></th>
            <th>PVIN model for time series simulations</th>
        </tr>
    </table>
</p>

### MATLAB

<p>
  <table style="width:100%">
        <tr>
            <th><code>PVIN_STEP_Simu.m</code></th>
            <th>Example code to generate step current-stimulated voltage response of PVINs, as in *Figs. 2,3 & 5* in *Ma et al. 2023*</th>
        </tr>
        <tr>
            <th><code>PVIN_Syn_Simu.m</code></th>
            <th>Example code to generate Poisson-distributed synaptic current-stimulated voltage response of PVINs, as in *Fig. 9* in *Ma et al. 2023*</th>
        </tr>    
        <tr>
            <th><code>PVIN_TwoParBifur.m</code></th>
            <th>Example code to reproduce the two-parameter bifurcation diagram shown in *Fig. 8* in *Ma et al. 2023*</th>
        </tr>      
        <tr>
            <th><code>NeuralCircuit.mlx</code></th>
            <th>Example code to simulate the spinal dorsal neural circuit model, as in *Fig. 9* in *Ma et al. 2023*</th>
        </tr>
        <tr>
            <th><code>PVIN_HH.m</code></th>
            <th>ODEs of PVIN model</th>
        </tr>
        <tr>
            <th><code>ePKCmodel.m<./code></th>
            <th>ODEs of ePKC&#947 interneuron model</th>
        </tr>      
        <tr>
            <th><code>PVIN_Cai.m</code></th>
            <th>Function: MATCONT equations of PVIN model</th>
        </tr>
        <tr>
            <th><code>PVIN_Cai.mat</code></th>
            <th>Function: MATCONT variables of PVIN model</th>
        </tr>
        <tr>
            <th><code>runHHmodel_STEP.m</code></th>
            <th>Function: Step current stimulation protocol of PVIN model</th>
        </tr>
        <tr>
            <th><code>runHHmodel_AbetaPoisson.m</code></th>
            <th>Function: Synaptic current stimulation protocol of PVIN model</th>
        </tr>
        <tr>
            <th><code>mySynInput.m</code></th>
            <th>Function: Generate the time series of AMPA and NMDA conductances at customized firing rate</th>
        </tr>
        <tr>
            <th><code>poissonSpikeGen.m</code></th>
            <th>Function: Generate Poisson-distributed timing of spike trains</th>
        </tr>
        <tr>
            <th><code>factor_syn.m</code></th>
            <th>Function: Generate a factor to correct the exponential equation peaks at 1</th>
        </tr>
        <tr>
            <th><code>genSyn.m</code></th>
            <th>Function: Generate the time series of a customized synaptic conductance</th>
        </tr>
      </table>
  </p>

<h2>IV. Instruction</h2>
<h3>XPPAUT<h3>
<h4>&emsp;Simulate the voltage response of PVIN model to step current stimulation</h4>
<ol>
    <li>Open <code>PVIN_2023.ode</code> with XPPAUT</li>
    <li>Parameter specification: injection current value <br>
        click 'Parameters', input 'iapp', input the value (e.g., 200), press 'Enter'
    </li>  
     <li>Parameter specification: calcium buffer concentration <br>
        click 'Parameters', input 'bt', input the value (e.g., 90), press 'Enter'
    </li>
    <li>Run the simulation  <br>
        click 'Initialconds' + '(G)o', and the result is like below <br> 
        <img src="./Figures/example_xpp.png" alt="PVIN model step current (xpp)" style="width:560px;height:400px;">
    </li>
</ol>
<h3>Matlab<h3>
<h4>&emsp;Simulate the voltage response of PVIN model to step current stimulation</h4>
<ol>
    <li>Run <code>PVIN_STEP_Simu.m</code></li>
    <li>In the command window, input the condition: 1 - naive PVIN | 2 - CCI PVIN </li>
    <li>In the command window, input the value of applied step current value (pA) </li>
    <li>A figure like below appears. The top panel shows voltage simulation and the bottom one is the applied current.</li>
</ol>
 <p><img src="./Figures/example_step.png" alt="PVIN model step current stimulation" style="width:560px;height:400px;"></p>

<h4>&emsp;Simulate the voltage response of PVINs to synaptic current stimulation</h4>
<ol>
    <li>Run <code>PVIN_Syn_Simu.m</code></li>
    <li>In the command window, input the condition: 1 - naive PVIN | 2 - CCI PVIN </li>
    <li>In the command window, input PVIN presynaptic current firing rate (Hz)  </li>
    <li>A figure like below appears. The top panel shows voltage simulation and the bottom one is the presyaptic conductances.</li>
</ol>
<p><img src="./Figures/example_syn.png" alt="PVIN model synaptic current stimulation" style="width:560px;height:400px;"></p>

<h4>&emsp;Reproduce the two-parameter bifurcation diagram shown in Fig. 8</h4>
<ol>
    <li>Run <code>PVIN_TwoParBifur.m</code></li>
    <li>After a few seconds, the two-parameter bifurcation diagram of PVIN model appears</li>
</ol>
<p><img src="./Figures/example_twopar.png" alt="PVIN model two-parameter bifurcation analysis" style="width:600px;height:600px;"></p>

<h4>&emsp;Simulate the PVIN embedded neural circuit model</h4>
<p>Open and run <Code>NeuralCircuit.mlx</Code>. This is a livescript presenting the simulation steps of the neural circuit model associated with Fig. 9 in Ma et al. 2023. The steps include: <br><p>
<ol>
    <li>Create Poisson-distributed synaptic conductances representing the presynaptic input from A&#946 fibers </li>
    <img src="./Figures/example_nc_abeta.png" alt="neural circuit model: Abeta fibersr" style="width:600px;height:a00px;">
    <li>Simulate the PVIN model using synaptic currents from A&#946 fibers</li>
    <img src="./Figures/example_nc_pvin.png" alt="neural circuit model: PVIN" style="width:600px;height:100px;">   
    <li>Simulate the PKC&#947 model using synaptic currernts from A&#946 fibers and PVIN model</li>
    <img src="./Figures/example_nc_pkcg.png" alt="neural circuit model: PKCgammaIN" style="width:600px;height:100px;">
</ol>

