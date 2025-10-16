a.


$\pi_{\text{name}} \left( \sigma_{\text{publisher} = \text{MCdraw-Hill}} (\text{books}) \bowtie_{\text{isbn}} \text{borrowed} \bowtie_{\text{memb_no}} \text{member} \right)$


b.


$\pi_{\text{name}} \left( \text{member} \where \text{memb_no} \in \left( \left( \sigma_{\text{publisher} = \text{MCdraw-Hill}} (\text{books}) \right) \div \left( \pi_{\text{isbn}} (\text{borrowed}) \right) \right) \right)$


c.


$\pi_{\text{name}, \text{memb_no}} \left( \text{member} \where \exists \geq 5 \left( \sigma_{\text{publisher} = \text{MCdraw-Hill}} (\text{books}) \bowtie \text{borrowed} \where \text{memb_no} = \text{memb_no} \right) \right)$


d.


$\pi_{\text{name}, \text{memb_no}, \text{publisher}} \left( \text{books} \bowtie \text{borrowed} \bowtie \text{member} \where \text{publisher} \land \text{COUNT}_\text{DISTINCT}(\text{isbn}) > 5 \right)$


e.


$\frac{\text{COUNT}(\text{borrowed})}{\text{COUNT}(\text{DISTINCT memb_no})} \text{ 在 } \text{borrowed} \text{ 中}$

渲染一下公式