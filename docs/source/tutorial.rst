===============
Tardis Tutorial
===============

lorem ipsum::

	$ some bash command here

lorem ipsum **bold lorem ipsum**.

-------------------------------
Tardis More Tutorial
-------------------------------

.. figure:: ../images/sample.png
   :alt: MoreFigure

   More Figure


lorem ipsum

1)	**bold**: 

	lorem ipsum

2)	**bold lorem ipsum**:

	lorem ipsum


----------------------
Tardis with some table
----------------------

lorem ipsum *italic lorem ipsum*

+-------------------------+---------------+-------------------+--------------+----------------------+
| Description             | Parameter     | Optional/Required | Precondition | Default Parameter    |
+=========================+===============+===================+==============+======================+
| Arg1                    | `-a`_         |  Optional         | None         | None                 |
+-------------------------+---------------+-------------------+--------------+----------------------+
| Arg2                    | `-b`_         |  Required         | `-a`_        | None ("path/to/file")|
+-------------------------+---------------+-------------------+--------------+----------------------+
| Args3                   | `-c`_         |  Required         | `-b`_        | `-a`_                |
|                         +---------------+                   |              |                      |
|                         | `-d`_         |                   |              |                      |
+-------------------------+---------------+-------------------+--------------+----------------------+
| Arg4                    | `-e`_         |  Required         | `-c`_        | None ("path/to/file")|
+-------------------------+---------------+-------------------+--------------+----------------------+
| Args5                   | `-f`_         |  Optional         | `-c`_        | None                 |
|                         +---------------+-------------------+--------------+----------------------+
|                         | `-g`_         |  Optional         | `-c`_        | None                 |
|                         +---------------+-------------------+--------------+----------------------+
|                         | `-h`_         |  Required         | `-d`_        | None ("path/to/file")|
|                         +---------------+-------------------+--------------+----------------------+
|                         | `-i`_         |  Required         | `-c`_        | `-c`_                |
|                         +---------------+                   |              |                      |
|                         | `-j`_         |                   |              |                      |
|                         +---------------+                   |              |                      |
|                         | `-k1`_        |                   |              |                      |
+-------------------------+---------------+-------------------+--------------+----------------------+


------------------------------------
Command-Line Parameters Descriptions
------------------------------------

Command line descriptions as shown in `Tardis with some table`_

-a
^^

lorem ipsum

-b
^^

lorem ipsum

-c
^^

lorem ipsum

-d
^^

lorem ipsum

-e
^^

lorem ipsum

-f
^^

lorem ipsum

-g
^^

lorem ipsum :option:`-c` is set. lorem ipsum

-h
^^

**Required** if :option:`-a` is set. lorem ipsum. See also `-k1`_.

-i
^^

lorem ipsum

-j
^^

lorem ipsum

-k1
^^^

lorem ipsum

	