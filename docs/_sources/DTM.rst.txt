Two methods for analyzing the communities on protein surface in DTM
===================================================================

Retrospective method
--------------------

When the approximate binding surface is known, it can be used to speed up the retrospective analysis.

1.	Identify the community overlapping with the binder in the binder complex.

2.	Run AlphaSpace on all snapshots and generate community based on each snapshot. And select the community by its
contact to the known binder, which was superimposed to each snapshot.

3. This gives the community of interests and designate it as the D-community, subsequent D-pocket analysis are done
only on the member of these communities.

This method, though simplistic, requires the information of the binding region, which is not always available. But for cases such as optimization of lead compound it provides an intuitive way represent the area of interests on the surface.


Perspective method
------------------

When the binding region is unknown and we would like to find the most targetable region on the protein surface.

1.	Simulation are carried out for apo protein to sample surface conformations.

2.	AlphaSpace are run on all snapshots and pockets are identified.

3.	All pockets are clustered based on lining atom Jaccard index, and each cluster is defined as D-pocket in DTM.

4.	The D-pockets are evaluated based on their per pocket space distribution, and based on a cutoff set either by
percentage of all D-pockets or hard cutoff, some are marked as Core-D-pocket.

5.	The list of Core-D-pocket gives a larger set of pockets in several snapshots, which will not all satisfy the
space requirements to be identified as Core-Pocket by itself, but are marked as Core-pocket anyway. This step is carried out to alleviate the fluctuations caused by a hard cutoff, so outlier pockets sharing similar trait would remain.

6.	From the set of identified Core-pocket, going back to their parent snapshot, we can identify the
connected auxiliary pockets by common lining atom count. This step expands around the Core-D-pocket and include the adjacent smaller pockets into the community. These pockets are marked Aux-pocket.

7.	Going back to the D-pocket list and examine the Non-Core-D-pockets. If a certain percentage of
pockets in a given Non-Core-D-pockets are marked Aux-pocket, that is they are connected to the Core-pockets in the corresponding snapshot, this D-pocket is marked as Aux-D-pocket. Unlike the previous method, here the Aux-D-pockets are defined so the communities they form with the Core-D-pockets would be more stable.