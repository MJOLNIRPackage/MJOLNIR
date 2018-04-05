Geometry
========

Description of the calculations converting the virtual positions of the Instrument object into energy transfer and A4 angles.

In this derivation it is assumed that the analysers are flat and have no focusing in the x-y direction when referring to the sketch in figure \ref{fig:SketchOfTransformationMatrixCalculation} but focusing in the x-z direction is allowed. In order to find the outgoing energy :math:`E_f` as well as the :math:`\vec{k}_f` wave vector, the position of the detector as seen from the sample is calculated. The energy is found from the :math:`A6` angle of the analyser. As the analyser can only change the momentum of the neutron in the plane of scattering, x-z plane, the path of the neutron in the x-y plane is to be a straight line as seen in the left of figure \ref{fig:SketchOfTransformationMatrixCalculation}. 

.. \begin{figure}[H]
.. \includestandalone[width=\linewidth]{Figures/SoftwareDocumentation/TransformationMatrices}
.. \caption{Sketch of detector positions relative to analysers and sample used in calculation of position in reciprocal space. \textbf{Inset} Overview of neutron path from sample to detector. \textbf{left}: 1:2 scaled path drawn on manifold following neutron path with $\vec{x}\prime=\cos{2\theta_A}\vec{x}+\sin{2\theta_A}\vec{z}$, such that the triangle relation is shown. \textbf{Left}: Neutron path from sample to detector as seen from above.}\label{fig:SketchOfTransformationMatrixCalculation}
.. \end{figure}

The scattering position on the analyser can be found from noting that the ratios,

.. math::

   \frac{\Delta X_D}{L_A+L_D}=\frac{\Delta X_A}{L_A}\qquad\Rightarrow\qquad\Delta X_A = \frac{\Delta X_D}{\frac{L_D}{L_A}+1}
 
are equal. From :math:`\Delta X_A` the actual scattering angle at the sample is found as

.. math::

   \cos{\tilde{\theta}_A}=&\frac{\tilde{\vec{L}}_A\cdot\tilde{\vec{L}}_D}{\left|{\tilde{\vec{L}}_A}\right|\left|{\tilde{\vec{L}}_D}\right|},\\
   \tilde{\vec{L}}_A=\begin{pmatrix}
   L_A\\ \Delta X_A \\ 0
   \end{pmatrix},&\qquad \tilde{\vec{L}}_D=\begin{pmatrix}
   L_D\cos{2\theta_A} \\ \Delta X_D\\L_D\sin{2\theta_A}
   \end{pmatrix}


Remembering \Linker{from section \ref{sec:Scattering}} that

.. math::

   \sin{\theta_A} = \frac{\pi}{k_fd} = \frac{\hbar\pi}{d \sqrt{2mE}},

where :math:`d` is the inter planar distance of the analyser. Thus, knowing the position of each detector pixel relative to the detector center one can find the correct scattering angle and final energy. 

.. \begin{wrapfigure}[20]{hr}{0.5\textwidth}\vspace{-0.4cm} 
.. %\begin{minipage}{\linewidth}
.. \centering
.. \includegraphics[trim={0cm 0cm 1.3cm 1.3cm},clip,width=0.95\linewidth]{Figures/SoftwareDocumentation/EfDetBank0}\\ %\label{fig:Ef}
.. \includegraphics[trim={0cm 0cm 1.3cm 1.3cm},clip,width=0.95\linewidth]{Figures/SoftwareDocumentation/ScatteringVectorDetBank0}
.. \caption{\textbf{Top}: Energy coverage and \textbf{bottom}: Scattering vector for all detector pixels for all three detectors in detector bank 0 and time bin 680 when calibrated using Vanadium measurements.}
.. \label{fig:EfScatteringVector}
.. \end{wrapfigure}

Next up is the initial energy which is found from simple distance-velocity calculations. The time taken from creation of neutron at the source to the scattering at the sample is given simply by distance divided by travel speed, or in energy

.. math::

   T_S=\frac{L_S}{\sqrt{\frac{2E_i}{m}}}

where the non-relativistic relation :math:`E=\frac{1}{2}mv^2` has been used and the distance travelled is denoted :math:`L_S`. The travelling time from sample to detector is

.. math::

   L_{SD} = \left|{\tilde{\vec{L}}_A}\right|+\left|{\tilde{\vec{L}}_D}\right|,

