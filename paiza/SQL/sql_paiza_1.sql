select h.id,h.name as name ,e.name as element,g.name as grade
from Hell as h
join Element as e
on h.element_id = e.id
join Grade as g
on g.id = h.grade_id
join ElementCompatibility as ec
on h.element_id = ec.weakness_element_id
join Hell as h2
on h2.name = "Graffiacane" and h2.element_id = ec.element_id
where g.name = "Boss" 