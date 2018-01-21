# Wireless-Network-Packet-Scheduling-Algorithms

Real time internet packets are deadline constrained. The packets are weighted to characterize relative importance over one another. The goal is to find an algorithm that maximizes weighted packet transmission before deadline expiry or satisfies quality of service for each link. This work develops discrete time simulation programs for various packet scheduling algorithms.

The model under consideration is interference constrained communication links, could be collocated or ad-hoc. On the other hand, packet arrival process is a Bernoulli random variable, independent from link to link. Please see [the presentation](https://drive.google.com/file/d/1cS8rGfJcWyW1awyca1oilJBFBPH_8e7j/view?usp=sharing) for a comprehensive idea about the project.

The following algorithms are studied-

# EDF: Earliest Deadline First
The algorithm prioritizes links having least remaining time packets. If two links have packets with the same least remaining time, then tie is broken giving priority to higher weight links.

# LDF: Largest Deficit First
Please see [the paper](http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6923479) for definition of deficit evoloution process. The algorithm schedules packets from the link that has highest deficit. Any possible tie is broken according to class weights.

# CMTO: Continuous Minloss Throughput Optimal
The algorithm looks at the existing packets waiting in the queue. It then finds the schedule that minimizes the loss over the queued packets. Please see [the paper](https://link.springer.com/article/10.1023%2FA%3A1015890105072) for detailed description.

# 3/4 Competitive Ratio Algorithm
This algorithm is an improvement over the previous ones. It achieves 75% of the optimal weighted packet transmission. Please see [the paper](http://www.columbia.edu/~js1353/pubs/jlss.pdf) for details.
