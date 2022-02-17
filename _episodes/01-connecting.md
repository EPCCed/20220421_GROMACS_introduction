---
title: "Connecting to ARCHER2 and transferring data"
teaching: 30
exercises: 15
questions:
- "How can I access ARCHER2 interactively and transfer data?"
objectives:
- "Understand how to connect to ARCHER2."
- "Know how to transfer data onto and off of ARCHER2 efficiently."
keypoints:
- "ARCHER2's login address is `login.archer2.ac.uk`."
- "The password policy for ARCHER2 is well documented."
- "There are a number of ways to transfer data to/from ARCHER2."
---

## Connecting using SSH

The login address for the full 23-cabinet ARCHER2 system is

```
login.archer2.ac.uk
```
{: .language-bash}

Access to ARCHER2 is via SSH using **both** a password and a passphrase-protected SSH key pair.

## Passwords and password policy

When you first get an ARCHER2 account, you will get a single-use password from the 
SAFE which you will be asked to change to a password of your choice. Your chosen 
password must have the required complexity as specified in the
[ARCHER2 Password Policy](https://www.archer2.ac.uk/about/policies/passwords_usernames.html).

The password policy has been chosen to allow users to use both complex, shorter passwords and
long, but comparatively simple passwords. For example, passwords in the style of both
`LA10!£lsty` and `horsebatterystaple` would be supported.

> ## Picking a good password
> Which of these passwords would be a good, valid choice according to the ARCHER2 Password
> Policy?
> 
> 1. `mypassword`
> 2. `rainbowllamajumping`
> 3. `A!94ufskl$?`
> 4. `horsebatterystaple`
> 
> > ## Solution
> >
> > 1. **No** This would not be accepted or a good choice as it is too short and is made up of obvious words
> > 2. **Yes** This would be a good choice as it is long enough and easy to remember
> > 3. **Yes** This would be accepted but may be difficult to remember and type (though you could use a password manager to store it)
> > 4. **No** While this meets the criteria, it is a well known example from a popular web comic and so would not be accepted
> >
> {: .solution}
{: .challenge}

## SSH keys

As well as password access, users are required to add the public part of an SSH key pair to access ARCHER2.
The public part of the key pair is associated with your account using the SAFE web interface.
See the ARCHER2 User and Best Practice Guide for information on how to create SSH key pairs
and associate them with your account:

* [Connecting to ARCHER2](https://docs.archer2.ac.uk/user-guide/connecting/)

## Data transfer services: scp, rsync, Globus Online

ARCHER2 supports several data transfer mechanisms. The one you choose depends
on the amount and structure of the data you want to transfer and where you want to transfer
the data to. The three main options are:

* `scp`: The standard way to transfer small to medium amounts of data (in the order of MBs or GBs) off ARCHER2 to any other location
* `rsync`: Used if you need to keep small to medium datasets (in the order of MBs or GBs) synchronised between two different locations
* *Globus Online*: Used to transfer large amounts of data (in the order of TBs) to other sites which are Globus Online enabled

More information on data transfer mechanisms can be found in the ARCHER2 User and Best Practice Guide:

* [Data management and transfer](https://docs.archer2.ac.uk/user-guide/data/)

## Data transfer best practice

Information on how to transfer data efficiently using the above methods is available in the ARCHER2 documentation. Nevertheless, here are the main points you should consider:

* **Not all data are created equal, understand your data.** Know what data you have. What is your
  critical data that needs to be copied to a secure location? Which data do you need in a different
  location to analyse? Which data would it be easier to regenerate rather than transfer? You should
  create a brief data management plan laying this out as this will allow you to understand which
  tools to use and when.
* **Minimise the data you are transferring.** Transferring large amounts of data is costly in both
  researcher time and actual time. Make sure you are only transferring the data you need to transfer.
* **Minimise the number of files you are transferring.** Each individual file has a static overhead in
  data transfers so it is efficient to bundle multiple files together into a single large
  archive file for transfer.
* **Does compression help or hinder?** Many tools have the option to use compression (e.g. `rsync`,
  `tar`, `zip`) and generally encourage you to use them to reduce data volumes. However, in some cases,
  the time spent compressing the data can take longer than actually transferring the uncompressed
  data; particularly when transferring data between two locations that both have large data transfer
  bandwidth available.
* **Be aware of encryption overheads.** When transferring data using `scp` (and `rsync` over `scp`)
  your data will be encrypted introducing a static overhead per file. This issue can be minimised by
  reducing the number files to be transferred by creating archives. You can also change the encryption
  algorithm to one that involves minimal encryption. The fastest performing cipher that is commonly 
  available in SSH at the moment is generally `aes128-ctr` as most common processors provide a
  hardware implementation.

> ## Creating an uncompressed zip archive and verifying the contents
> Using the documentation above, find the command you would use to create an uncompressed zip archive
> file of all data within a directory called `large_data_output/`. What command would you use to verify
> that the archive file you have created is not corrupt so you can safely delete the original data?
> > ## Solution
> > You use the `zip` command to archive the data. The `-r` option is used to perform the operation
> > recursively on a directory and the `-0` option is used to specify the archive should be uncompressed:
> > ```
> > auser@ln01:~> zip -0r large_data_output.zip large_data_output/
> > ```
> > {: .language-bash}
> > To verify the archive is valid, you would use the `zip` command again, this time with the `-T` 
> > option:
> > ```
> > auser@ln01:~> zip -T large_data_output.zip
> > ```
> > {: .language-bash}
> {: .solution}
{: .challenge}

{% include links.md %}

