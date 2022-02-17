---
title: "Using SAFE"
teaching: 30
exercises: 0
questions:
- "How can I manage my ARCHER2 user account?"
- "How can I gain access to licensed software on ARCHER2?"
- "How can I check what resources I am using and look at historical usage?"
objectives:
- "Know what SAFE can be used for and how to use it."
- "Understand what reports you can generate on your use of ARCHER2."
keypoints:
- "SAFE allows you to manage accounts and resources online."
---

SAFE is a web-based tool developed by EPCC that allows you to:

* Manage your ARCHER2 account
* Record research outputs linked to your use of ARCHER2
* Look at your historical use of ARCHER2 service
* (If you are a project leader/manager) Manage your project's resources on ARCHER2

You should already have used SAFE to apply for your ARCHER2 user account for this course. In
this section we will look at how you use SAFE for common tasks associated with your use of
ARCHER2.

More information on using SAFE can be found in
[the SAFE User Guide](https://epcced.github.io/safe-docs/)

## Account management

You will have at least two accounts associated with your use of ARCHER2: an account on the
SAFE web tool and ARCHER2 login (sometimes called *machine*) accounts.

> ## Multiple login accounts
> You can have multiple ARCHER2 login accounts associated with a single SAFE account. Many
> users will also have access to the EPSRC Tier-2 HPC resources and login accounts for those
> systems will often also be linked to our SAFE account. So you can actually have multiple
> login accounts on multiple services associated with a single SAFE account!
{: .callout}

We will not document the full use of SAFE here but will mention a few of the most useful
options with managing your different accounts.

### Your SAFE account

Your SAFE account is used to allow you to manage your login accounts and resources, set
contact details and preferences, record research outputs associated with your use of the
service and report on your historical use of the service. We cover these topics in more
detail below but first we will look at how to use your institutional single sign on (SSO)
credentials to sign into SAFE and how to update your e-mail notification preferences.

**Using your institutional SSO to log in to SAFE**. You should first log into SAFE using
the e-mail and password combination that you would normally use to log in, then use the
*Your details - Register institutional ID* menu option. This should take you to a dialog
box where you can select your home institution. Once you select your institution, you will
be taken to your institution SSO page (if you are not already signed in) and this will then
be linked to your SAFE account. You will now be able to use the *Login with institutional
credentials* button on the SAFE login page to login to the SAFE using your institutional SSO
ID.

**Changing your e-mail notification preferences**. There are a number of different mailing
lists associated with the ARCHER2 service and you can subscribe to any combination of these.
You can change your preferences in SAFE by use the *Your details - Update email settings*
menu option. You can then select the mailing lists you wish to subscribe to.

### Your ARCHER2 login account(s)

There can be multiple login accounts associated with a single SAFE account. Each of these
accounts will be listed in the *Login accounts* menu. Selecting a particular login account
in the menu will give you access to a page with more options for managing those accounts.
We will not cover all options here, but briefly mention a few key options.

**Viewing resources available**. When you select a login account you will see a summary
of the systems it has access to, the resources (both compute allocations and storage quotas)
available and used, and lists of licensed software that the account has access to.

**Picking up your initial login password**. When your login account is first created you
will pick up your initial, single-use password from the *View login password* button on
the appropriate login account page. When you use this password to log in to ARCHER2, you
will be immediately asked to change this.

> ## SAFE is not kept up to date with your actual password
> Note that the login account password stored in SAFE is not kept up to date with your
> password on the ARCHER2 machine once you have changed it. If you forget your password,
> you will need to reset it as described below.
{: .callout}

**Resetting your login password**. If you forget your login password or need to reset 
it for any other reason, you can use the *Request Password Reset* button to ask for
a new single-use password to log into the service. You will be notified by e-mail once
the new password is available for you to pick up and you use it exactly the same way
as your initial password.

**Adding SSH keys to your account**. To log into ARCHER2 you need to use both a password
and an key pair SSH protected by a strong passphrase. Once you have generated your 
SSH key pair you use SAFE to associate the public part of the key pair with your account.
Use the *Add credential* button and then select the *SSH key* option to associate the
public part of the key pair with your account.

**Requesting access to licensed software**. If you wish to access particular licensed
software on ARCHER2 you can request access by using the *Request Access to Package*
button on the Login account page. This will take you to a form where you can provide 
you licence details that allow you to access the software. The ARCHER2 Service Desk
will verify these details (sometimes by contacting the software developers or licence
contact) and then enable access.

## Recording research outputs

In order to make the case to UK Government for continued support of supercomputing and
HPC we need to collect research outputs from all the projects. The SAFE provides a 
mechanism for uploading DOI's and automatically gathering the output details based on
the DOI's provide. You can also get reports out of SAFE to assist with submission to,
for example, ResearchFish and for project reporting. You add DOI's to your choice of 
your projects in SAFE by using the *Your details - Publications* menu item.

## Reporting on your use of the service

SAFE also provides a way to investigate your historical use of the service via the
Reporting functionality. This is accessed via the *Service information - Reports*
menu item. There are various reports available that provide details on single 
jobs up to overviews of usage over particular periods.

Output from all reports is available in different formats (HTML, PDF and CSV).
Usually, you will run them using the *Preview* button to display the output in the
web browser and provide the interface to change report parameters so you can 
explore your use of the system.

<!-- TODO: Add exercise on running SAFE report on usage once these are available.
Will SAFE get the data from jobs run on the course soon enough to be able to use
this as an example? If not, use a generic report. -->

{% include links.md %}

