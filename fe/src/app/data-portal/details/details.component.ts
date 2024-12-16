import {Component, OnInit} from '@angular/core';
import {ActivatedRoute} from "@angular/router";
import {ElasticService} from "../../services/elastic.service";
import {MatListModule} from "@angular/material/list";

@Component({
  selector: 'app-details',
  standalone: true,
  imports: [MatListModule],
  templateUrl: './details.component.html',
  styleUrl: './details.component.css'
})
export class DetailsComponent implements OnInit {
  data: any;
  constructor(
    private activatedRoute: ActivatedRoute,
    private _elasticService: ElasticService
  ) {}

  ngOnInit() {
    this._elasticService.getRecordDetails(this.activatedRoute.snapshot.paramMap.get('urn')).subscribe(
      data => {
        this.data = data.results[0];
      }
    );
  }
}
