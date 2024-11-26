# Combinatoric_background_e1039
Mixing method explained in (https://arxiv.org/abs/2302.04152) was implemented by Abinash Pun for e906 experiment. A modified version that can extract data from e1039 DST files and perform new SQVertexing is included here.

## Key modifications 
Two major changes were made to make the scripts compatible with the e1039 DST files.
- E906 had SRawEvent and SRecEvent objects in the DST file. But SQEvent and SRecEvent for e1039. Created the SRawEvent objects using SQEvent and SQHitVector objects available in e1039 DST file
- E906 used old `VertexFit` algorithm to perform the vertexing of the tracks. But e1039 data is reconstructed using new `SQVeretexing` algorithm. A custom version (outside the Fun4All framework) was created to introduce the new vertxign to this project.
