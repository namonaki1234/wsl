select h.id,h.name as name ,e.name as element,g.name as grade
from Hell as h
join Element as e
on h.element_id = e.id
join Grade as g
on g.id = h.grade_id
where e.name = "Air" and g.name = "Boss"