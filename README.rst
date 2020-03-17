
.. image:: https://github.com/3DGenomes/tadbit/raw/master/doc/source/pictures/TADbit_logo.png
   :height: 50
   :width: 240



+-------------------------------------+---------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+---------------------------------------------------------------+
|                                     | .. image:: https://travis-ci.org/3DGenomes/TADbit.png?branch=master       | .. image:: https://coveralls.io/repos/github/3DGenomes/TADbit/badge.svg?branch=master       | .. image:: https://img.shields.io/badge/license-GPL-green.svg |
| Current version: v0.4.97            |   :target: https://travis-ci.org/3DGenomes/TADbit                         |   :target: https://coveralls.io/github/3DGenomes/TADbit?branch=master                       |                                                               |
|                                     |                                                                           |                                                                                             |                                                               |
+-------------------------------------+---------------------------------------------------------------------------+---------------------------------------------------------------------------------------------+---------------------------------------------------------------+


TADbit is a complete Python library to deal with all steps to analyze,
model and explore 3C-based data. With TADbit the user can map FASTsQ
files to obtain raw interaction binned matrices (Hi-C like matrices),
normalize and correct interaction matrices, identify and compare the
Topologically Associating Domains (TADs), build 3D models
from the interaction matrices, and finally, extract structural
properties from the models. TADbit is complemented by TADkit for
visualizing 3D models.

Hi-C experiments generate genomic interaction between loci located in
the same or in different chromosomes. TADbit is built around the
concept of a chromosome, and uses it as a central item to store and
compare different Hi-C experiments. The library has been designed to
be used by researchers with no expertise in computer
science. All-in-one scripts provided in TADbit allow to run the full
analysis using one single command line; advanced users may produce
their own programs using TADbit as a complementary library.


Contributors
************

TADbit is currently developed at the  `MarciusLab <http://www.marciuslab.org>`_ with the contributions of François Serra, David Castillo, Marco Di Stefano, Irene Farabella, Mike Goodstadt and many other members of our Lab

Documentation
*************

* `TADbit documentation <http://3dgenomes.github.io/TADbit/>`_
* `Summary of TADbit classes and functions <https://github.com/3DGenomes/tadbit/blob/master/doc/summary.rst>`_

Feedback
--------
If you have any question remaining, we would be happy to answer informally:

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/3DGenomes/tadbit
   :target: https://gitter.im/3DGenomes/tadbit?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

Frequently asked questions
--------------------------

Check the label `FAQ <https://github.com/3DGenomes/TADbit/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3AFAQ+>`_ in TADbit issues.

If your question is still unanswered feel free to open a new issue.

Docker/Singularity Containers
-----------------------------

Recipe files (`Dockerfile <https://docs.docker.com/engine/reference/builder/>`_ and
`Singularity recipe <https://www.sylabs.io/guides/2.6/user-guide/quick_start.html#build-images-from-scratch>`_) to generate containers are
available in the containers folder.


- **Docker**

Build the image using the :code:`Dockerfile` from inside an empty folder with :code:`docker build -t tadbit .` (~20 minutes)

Once built, run it as :code:`docker run tadbit tadbit map -h`

This image contains all dependencies for TADbit and also `jupyter <http://jupyter.org/>`_ .

To run a notebook from inside the docker container run :code:`tadbit` docker image as::

    docker run -it -p 8888:8888 -v /LOCAL_PATH:/mnt tadbit

:code:`LOCAL_PATH` *would be for example a local folder with data*
*(e.g. FASTQs or reference genomes). And* :code:`/mnt` *a directory*
*inside the Docker container where the* :code:`LOCAL_PATH` *would be mounted.*

From inside docker run::

  jupyter notebook --ip 0.0.0.0 --allow-root --NotebookApp.token=''

And finally write the url :code:`http://localhost:8888` in your browser.

*Note: this can also be done in a single line and running in the background:*

::

  docker run -d -p 8888:8888 -v /LOCAL_PATH:/mnt tadbit jupyter notebook --ip 0.0.0.0 --allow-root --NotebookApp.token='' > /dev/null &

- **Singularity**

Build the image using the :code:`Singularity` from inside an empty folder with :code:`sudo singularity build tadbit.simg Singularity` (~20 minutes)

Once built, run it as :code:`singularity run tadbit.simg`

You can also install jupyter inside the Singularity by uncommenting the coresponding line in the recipe file.

Citation
********
Please, cite this article if you use TADbit.

Serra, F., Baù, D., Goodstadt, M., Castillo, D. Filion, G., & Marti-Renom, M.A. (2017).
**Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors.**
*PLOS Comp Bio* 13(7) e1005665. `doi:10.1371/journal.pcbi.1005665 <https://doi.org/10.1371/journal.pcbi.1005665>`_

Methods implemented in TADbit
-----------------------------
In addition to the general citation for the TADbit library, please cite these articles if you used TADbit for:

- Mapping and read filtering [Marco-Sola2012]_ [Imakaev2012]_ [Ay2015]_
- Hi-C normalization [Imakaev2012]_ [Rao2014]_
- A/B compartment calling [Lieberman-Aiden2009]_
- Model assessement [Trussart2015]_
- Chromatin 3D Model Building [BaùMarti-Renom2012]_

Applications
------------
TADbit has been previously used for modeling genomes and genomic domains. Here is the list of published articles:

- Alpha-globin domain [Baù2011]_
- *Caulobacter crescentus* genome [Umbarger2011]_
- TADs as regulons [Le_Dily2014]_
- Yeast chromosome III [Belton2015]_
- *Mycoplasma pneumoniae* genome [Trussart2017]_
- Comparison with super-res imaging [Cattoni2017]_
- Analysis of time-series Hi-C data during transdifferentiation [Stadhouders2018]_
- Analysis of cohesin subunits effect on the 3D genome [Kojic2018]_ [Cuadrado2019]_
- Analysis of MLL mutant on the 3DGenome [Mas2018]_
- Analysis of internal lamin on the 3DGenome [Pascual-Reguant2018]_
- Integration with IMGR of super-res data and Hi-C data [Nir2018]_
- 3D enhancer hubs in diabetes [Miguel-Escalada2019]_
- Modeling location of RNA in the nucleus [Morf2019]_

Other programs
--------------
TADbit uses other major software packages in biology. Here is the list of their articles:

- IMP Integrative Modeling Platform [Russel2011]_
- MCL Markov Cluster Algorithm [Enright2002]_

TADbit training
***************

Next editions
-------------

* To be announced.

Past editions
-------------

* April 10th to April 11th 2017: `MuG
  <http://www.multiscalegenomics.eu/MuGVRE/>`_ workshop: `Multi-scale study of 3D Chromatin structure
  <http://www.multiscalegenomics.eu/MuGVRE/multi-scale-study-of-3d-chromatin-structure/>`_ at the
  `European Bioinformatics Institute (EMBL-EBI) <https://www.embl.de/training/cco/>`_,
  Hinxton, Cambridge, (United Kingdom)
* April 3rd to April 7th 2017: `Chromosomal Conformation course
  <http://www.crg.eu/en/event/coursescrg-chromosomal-conformation-0>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Barcelona (Spain)
* October 10th to October 14th 2016: `3DAROC16 3C-based data analysis and 3D reconstruction of chromatin folding
  <http://gtpb.igc.gulbenkian.pt/bicourses/3DAROC16/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)
* September 28th to October 2nd 2015: `Chromosomal Conformation course
  <http://gtpb.igc.gulbenkian.pt/bicourses/2014/CSDM14/>`_ at the
  `CRG <http://www.crg.eu/en/content/training/>`_
  training programme Barcelona (Spain)
* November 25th to November 28th 2014: `CSDM 2014
  <http://gtpb.igc.gulbenkian.pt/bicourses/2014/CSDM14/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)
* September 6th 2014: `TADbit: Automated Analysis and
  Three-Dimensional Modeling of Genomic Domains
  <http://www.eccb14.org/program/tutorials/tadbit>`_ at `ECCB14
  <http://www.eccb14.org/>`_ Strasbourg (France)
* November 27th to November 29th 2013: `CSDM 2013
  <http://gtpb.igc.gulbenkian.pt/bicourses/2013/CSDM13/>`_ at the
  `GTPB <http://gtpb.igc.gulbenkian.pt/bicourses/index.html>`_
  training programme Oeiras (Portugal)


Bibliography
************


.. [Ay2015] Ay, F., Vu, T.H., Zeitz, M.J., Varoquaux, N., Carette, J.E., Vert, J.-P., Hoffman, A.R. and Noble, W.S. 2015. Identifying multi-locus chromatin contacts in human cells using tethered multiple 3C. BMC Genomics 16, p. 121.

.. [BaùMarti-Renom2012] Baù, D. and Marti-Renom, M.A. 2012. Genome structure determination via 3C-based data integration by the Integrative Modeling Platform. Methods 58(3), pp. 300–306.

.. [Baù2011] Baù, D., Sanyal, A., Lajoie, B.R., Capriotti, E., Byron, M., Lawrence, J.B., Dekker, J. and Marti-Renom, M.A. 2011. The three-dimensional folding of the α-globin gene domain reveals formation of chromatin globules. Nature Structural & Molecular Biology 18(1), pp. 107–114.

.. [Belton2015] Belton, J.-M., Lajoie, B.R., Audibert, S., Cantaloube, S., Lassadi, I., Goiffon, I., Baù, D., Marti-Renom, M.A., Bystricky, K. and Dekker, J. 2015. The conformation of yeast chromosome III is mating type dependent and controlled by the recombination enhancer. Cell reports 13(9), pp. 1855–1867.

.. [Cattoni2017] Cattoni, D.I.,  Cardozo-Gizz, A.M.,  Georgieva, M.,  Di Stefano, M.,  Valeri, A.,  Chamousset, D.,  Houbron, C.,  Dejardin, S.,  Fiche, J-B.,  Marti-Renom, M.A.,  Bantignies, F.,  Cavalli, G. and Nollmann, M. (2017) Single-cell absolute contact probability detection reveals that chromosomes are organized by modulated stochasticity. Nature Communications 8 pp 1753

.. [Cuadrado2019] Cuadrado, A.,  Giménez-Llorente, D.,  Kojic, A.,  Rodríguez-Corsino, M.,  Cuartero, Y.,  Martín-Serrano, G.,  Gómez-López, G.,  Marti-Renom, M.A. and Losada, A. (2019) Specific contributions of cohesin-SA1 and cohesin-SA2 to TADs and Polycomb domains in embryonic stem cells. Cell Reports, in press

.. [Enright2002] Enright, A. J., Van Dongen, S., & Ouzounis, C. A. (2002). An efficient algorithm for large-scale detection of protein families. Nucleic Acids Research, 30(7), 1575–1584.

.. [Imakaev2012] Imakaev, M., Fudenberg, G., McCord, R.P., Naumova, N., Goloborodko, A., Lajoie, B.R., Dekker, J. and Mirny, L.A. 2012. Iterative correction of Hi-C data reveals hallmarks of chromosome organization. Nature Methods 9(10), pp. 999–1003.

.. [Kojic2018] Kojic, A.,  Cuadrado, A.,  Koninck, A.M.,  Gomez-Lopez, G.,  Rodriguez-Corsino, M.,  Le Dily, F.,  Marti-Renom, M.A. and Losada, A. (2018) Distinct roles of cohesin-SA1 and cohesin-SA2 in 3D chromosome organization. Nature Structural and Molecular Biology 25 pp 496–504

.. [Le_Dily2014] Le Dily, F., Baù, D., Pohl, A., Vicent, G.P., Serra, F., Soronellas, D., Castellano, G., Wright, R.H.G., Ballare, C., Filion, G., Marti-Renom, M.A. and Beato, M. 2014. Distinct structural transitions of chromatin topological domains correlate with coordinated hormone-induced gene regulation. Genes & Development 28(19), pp. 2151–2162.

.. [Lieberman-Aiden2009] Lieberman-Aiden, E., van Berkum, N.L., Williams, L., Imakaev, M., Ragoczy, T., Telling, A., Amit, I., Lajoie, B.R., Sabo, P.J., Dorschner, M.O., Sandstrom, R., Bernstein, B., Bender, M.A., Groudine, M., Gnirke, A., Stamatoyannopoulos, J., Mirny, L.A., Lander, E.S. and Dekker, J. 2009. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. Science 326(5950), pp. 289–293.

.. [Marco-Sola2012] Marco-Sola, S., Sammeth, M., Guigo, R. and Ribeca, P. 2012. The GEM mapper: fast, accurate and versatile alignment by filtration. Nat Methods 9(12), pp. 1185-1188.

.. [Mas2018] Mas, G.,  Blanco, E.,  Ballaré, C.,  Sansó, M.,  Spill, Y.G.,  Hu, D.,  Aoi, Y.,  Le Dily, F.,  Shilatifard, A.,  Marti-Renom, M.A. and Di Croce, L. (2018) Promoter bivalency favors an open architecture of the stem cell genome. Nature Genetics 50 pp 1452–1462

.. [Miguel-Escalada2019] Miguel-Escalada, I.,  Bonàs-Guarch, S.,  Cebola, I.,  Ponsa-Cobas, J.,  Mendieta-Esteban, J. ,  Rolando, D.,  Javierre, B.M.,  Atla, G.,  Farabella, I.,  Morgan, C.C.,  García-Hurtado, J.,  Beucher, A.,  Morán, I.,  Pasquali, L.,  Ramos, M.,  Appel, E.V.R.,  Linneberg, L.,  Gjesing, A.P.,  Witte, D.R.,  Pedersen, O.,  Grarup, N.,  Ravassard, P.,  Mercader, J.M.,  Torrents, D.,  Piemonti, L.,   Berney, T.,  de Koning E.,  Kerr-Conte, J.,  Pattou, F.,  Hansen, T.,   Marti-Renom, M.A.,  Fraser, P. and Ferrer, J. (2019) Human pancreatic islet 3D chromatin architecture provides insights into the genetics of type 2 diabetes. Nature Genetics, in press

.. [Morf2019] Morf, J.,  Wingett, S.W.,  Farabella, I.,  Cairns, J.,   Furlan-Magaril, M.,  Jiménez-García, L.F.,  Liu, X.,  Craig, F.F.,  Walker, S.,  Segons-Pichon, A.,  Andrews, S.,  Marti-Renom, M.A. and Fraser, P. (2019) RNA proximity sequencing reveals properties of spatial transcriptome organization in the nucleus. Nature Biotechnology, in press

.. [Nir2018] Nir, G.,  Farabella, I.,  Pérez Estrada, C.,   Ebeling, C.G.,  Beliveau, B.J.,  Sasaki, H.M.,  Lee, S.H.,  Nguyen, S.C.,  McCole, R.B.,  Chattoraj, S.,  Erceg, J.,  Abed, J.A.,  Martins, N.M.C.,   Nguyen, H.Q.,  Hannan, M.A.,  Russell, S.,  Durand, N.C.,  Rao, S.S.P.,  Kishi, J.Y.,  Soler-Vila, P.,  Di Pierro, M.,  Onuchic, J.N.,  Callahan, S.,  Schreiner, J.,  Stuckey, J.,  Yin, P.,  Lieberman Aiden, E.,  Marti-Renom, M.A. and Wu, C.T. (2018) Walking along chromosomes with super-resolution imaging, contact maps, and integrative modeling. PLOS Genetics 14(12) pp e1007872

.. [Pascual-Reguant2018] Pascual-Reguant. L.,  Blanco, E.,  Galan, S.,  Le Dily, F.,  Cuartero, Y.,  Serra-Bardenys, G.,  di Carlo, V.,  Iturbide, A.,  Cebrià-Costa, J.P.,  Nonell, L.,  García de Herreros, A.,  Di Croce, L.,  Marti-Renom, M.A. and Peiró, S. (2018) Genome-wide mapping of lamin B1 reveals the existence of dynamic and functional euchromatin lamin B1 domains (eLADs) during epithelial-to-mesenchymal transition (EMT).Nature Communications 9(1) pp 3420

.. [Rao2014] Rao, S.S.P., Huntley, M.H., Durand, N.C., Stamenova, E.K., Bochkov, I.D., Robinson, J.T., Sanborn, A.L., Machol, I., Omer, A.D., Lander, E.S. and Aiden, E.L. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 159(7), pp. 1665–1680.

.. [Russel2011] Russel, D., Lasker, K., Webb, B., Velázquez-Muriel, J., Tjioe, E., Schneidman-Duhovny, D., et al. (2011). Putting the Pieces Together: Integrative Modeling Platform Software for Structure Determination of Macromolecular Assemblies. PLoS Biology, 10(1), e1001244.

.. [Stadhouders2018] Stadhouders, R.,  Vidal, E.,  Serra, F.,  Di Stefano, B.,  Le Dily, F.,  Quilez, J.,  Gomez, A.,  Collombet, S.,  Berenguer, C.,  Cuartero, Y.,  Hecht, J.,  Filion, G.,  Beato, M.,  Marti-Renom, M.A. and Graf, T. (2018) Transcription factors orchestrate dynamic interplay between genome topology and gene regulation during cell reprogramming. Nature Genetics 50 pp 238–249

.. [Trussart2015] Trussart, M., Serra, F., Baù, D., Junier, I., Serrano, L. and Marti-Renom, M.A. 2015. Assessing the limits of restraint-based 3D modeling of genomes and genomic domains. Nucleic Acids Research 43(7), pp. 3465–3477.

.. [Trussart2017] Trussart, M., Yus, E., Martinez, S., Baù, D., Tahara, Y.O., Pengo, T., Widjaja, M., Kretschmer, S., Swoger, J., Djordjevic, S., Turnbull, L., Whitchurch, C., Miyata, M., Marti-Renom, M.A., Lluch-Senar, M. and Serrano, L. 2017. Defined chromosome structure in the genome-reduced bacterium Mycoplasma pneumoniae. Nature Communications 8, p. 14665.

.. [Umbarger2011] Umbarger, M.A., Toro, E., Wright, M.A., Porreca, G.J., Baù, D., Hong, S.-H., Fero, M.J., Zhu, L.J., Marti-Renom, M.A., McAdams, H.H., Shapiro, L., Dekker, J. and Church, G.M. 2011. The three-dimensional architecture of a bacterial genome and its alteration by genetic perturbation. Molecular Cell 44(2), pp. 252–264.
