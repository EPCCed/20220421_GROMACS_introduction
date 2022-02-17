---
title: "Bootstrapping your use of ARCHER2"
teaching: 45
exercises: 90
questions:
- "How can I get started with using ARCHER?"
objectives:
- "Get help to get your work up and running on the ARCHER2 system"
keypoints:
- "Understand the next steps for you in using ARCHER2."
---

Now you know enough about ARCHER2 to explore how to use it for your work or to understand
what its potential benefits are you. You may also have ideas around where the 
barriers and difficulties may lie and have further questions on how you can 
start using and/or trying ARCHER2 for your work.

This session is designed to give you the opportunity to explore these questions and
issues. The instructors and helpers on the course will be on hand to answer your
questions and discuss next steps with you.

> ## Potential discussions
>
> Things you could discuss with the instructors and helpers could include:
>
> - Your computational workflow and where ARCHER2 could help
> - How to get access to ARCHER2 for your work
> - How to get help and support to get your work running using ARCHER2.
>   For example, software development, further training, access to local expertise
{: .callout}

## Options for this session

There are a number of different options for practical work during this session. The
challenges below include: exploring your own work; and an extended example using a parallel
HPC application. If you have something else you want to use the session for (e.g. to
discuss things with the instructors/helpers as described above) then please feel free
to do this. The idea of the session is to help you bootstrap your use of ARCHER2
and this will differ from individual to individual!

> ## Exploring your work using ARCHER2
>
> If you have a practical example of something from your area of work that you would like
> help with getting up and running on an ARCHER2 or exploring the performance of
> on an ARCHER2, this is great! Please feel free to discuss this with us and ask
> questions (both technical and non-technical).
{: .challenge}

> ## Exploring the performance of GROMACS
>
> [GROMACS](http://www.gromacs.org) is a world-leading biomolecular modelling package
> that is heavily used on HPC systems around the world. Choosing the best resources
> for GROMACS calculations is non-trivial as it depends on may factors, including:
>
> - The underlying hardware of the HPC system being used
> - The actual system being modelled by the GROMACS package
> - The balance of processes to threads used for the parallel calculation
>
> In this exercise, you should try and decide on a good choice of resources and settings
> on ARCHER2 for a typical biomolecular system. This will involve:
>
> - Downloading the [input file for GROMACS]({{site.github.repository_url}}/blob/gh-pages/files/ion_channel.tpr?raw=true)
> - Writing a job submission script to run GROMACS on ARCHER2 using the system documentation
> - Varying the number of nodes (from 1 to 16 nodes is a good starting point) used for the GROMACS job
>   and benchmarking the performance (in ns/day)
> - Using the results from this study to propose a good resource choice for this GROMACS calculation
>
> If you want to explore further than this initial task then there are a number of 
> different interesting ways to do this. For example:
> 
> - Vary the number of threads used per process
> - Reduce the number of cores used per node
> - Allow the calculation to use Symmetric Mutithreading (SMT)
>
> Please ask for more information on these options from a helper!
{: .challenge}


{% include links.md %}


